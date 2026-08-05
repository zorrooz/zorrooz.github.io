import OpenAI from 'openai'
import fs from 'fs/promises'
import fsSync from 'fs'
import path from 'path'
import crypto from 'crypto'
import { pathToFileURL } from 'url'
import { cacheDir, contentSrcDir } from '../dataConfig.ts'

let clientCache: OpenAI | null = null

/**
 * 懒加载 DeepSeek 配置（llmConfig.ts 被 gitignore，CI 无此文件）。
 * 与 tagMerger 相同的动态加载策略：模块可被 vite-ssg 打包（静态 import 会
 * 在 esbuild 阶段因缺文件失败），翻译只在本地/有配置时进行。
 */
async function loadLlmConfig(): Promise<{ url: string; apikey: string; model: string } | null> {
  const cfgPath = path.join(import.meta.dirname, '../../config/llmConfig.ts')
  if (!fsSync.existsSync(cfgPath)) return null
  try {
    const mod = (await import(pathToFileURL(cfgPath).href)) as {
      default?: { url?: string; apikey?: string; model?: string }
    }
    const cfg = mod.default || {}
    if (!cfg.url || !cfg.apikey || !cfg.model) return null
    return { url: cfg.url, apikey: cfg.apikey, model: cfg.model }
  } catch {
    return null
  }
}

async function getClient(): Promise<OpenAI> {
  if (clientCache) return clientCache
  const cfg = await loadLlmConfig()
  if (!cfg) throw new Error('未找到 LLM 配置：请在 src/config/llmConfig.ts 配置 { url, apikey, model }')
  clientCache = new OpenAI({ baseURL: cfg.url, apiKey: cfg.apikey })
  return clientCache
}

/**
 * 翻译文本内容
 * @param {string} text - 要翻译的文本
 * @param {string} fileType - 文件类型 ('md' 或 'yaml')
 * @returns {Promise<string>} 翻译后的文本
 */
async function translateText(text: string, fileType = 'md'): Promise<string> {
  try {
    const [client, cfg] = await Promise.all([getClient(), loadLlmConfig()])
    const systemPrompt =
      fileType === 'yaml'
        ? '你是一个专业翻译器。请将以下YAML文件内容翻译为英文，要求：\n1. 严格保持YAML格式结构（键名不翻译，只翻译值）\n2. 保持缩进、标点等格式不变\n3. 不添加任何解释、注释或额外内容\n4. 只输出翻译结果'
        : '你是一个专业翻译器。请将以下Markdown内容翻译为英文，要求：\n1. 严格保持原始格式（包括Markdown语法、代码块、换行、缩进等）\n2. 不添加任何解释、注释或额外内容\n3. 只输出翻译结果'

    const completion = await client.chat.completions.create({
      messages: [
        { role: 'system', content: systemPrompt },
        { role: 'user', content: text },
      ],
      model: cfg?.model ?? '',
      temperature: 0.3,
    })

    const usage = completion.usage
    if (usage) {
      console.log(
        `[TOKENS] prompt=${usage.prompt_tokens} completion=${usage.completion_tokens} total=${usage.total_tokens}`,
      )
    }

    return (completion.choices[0].message.content ?? '').trim()
  } catch (error) {
    console.error('翻译时发生错误:', (error as Error).message)
    throw error
  }
}

/**
 * 内容哈希（增量翻译判定依据）：比 mtime 可靠——
 * 文件被复制/同步/切换分支导致 mtime 变化时不会误触发重新翻译。
 */
function contentHash(content: string): string {
  return crypto.createHash('md5').update(content, 'utf-8').digest('hex')
}

/**
 * 状态键：相对 content-src 的路径（不依赖机器绝对路径），
 * 使 cache/.translate-state.json 可跨机器/跨分支复用。
 */
function toStateKey(sourcePath: string): string {
  return path.relative(contentSrcDir, sourcePath).split(path.sep).join('/')
}

interface TranslationState {
  [sourcePath: string]: string
}

function stateFilePath(_dir: string): string {
  // 翻译状态统一收敛到 cache/（数据分支），不再散落在各源目录
  return path.join(cacheDir, '.translate-state.json')
}

function loadState(dir: string): TranslationState {
  try {
    return JSON.parse(fsSync.readFileSync(stateFilePath(dir), 'utf-8')) as TranslationState
  } catch {
    return {}
  }
}

function saveState(dir: string, state: TranslationState): void {
  // 防御：空状态不写，避免覆盖已建立的翻译状态（否则下次会全量重翻）
  if (!state || Object.keys(state).length === 0) return
  try {
    fsSync.writeFileSync(stateFilePath(dir), JSON.stringify(state, null, 2), 'utf-8')
  } catch (e) {
    console.warn(
      `[Warn] 翻译状态写入失败（不影响本次翻译，但下次会重复翻译）: ${
        e instanceof Error ? e.message : e
      }`,
    )
  }
}

/**
 * 检查文件是否需要翻译（增量翻译）：目标不存在，或源文件内容哈希与上次翻译时不同。
 * @param {string} sourcePath - 源文件路径
 * @param {string} targetPath - 目标文件路径
 * @param {TranslationState} state - 翻译状态（源路径 → 内容哈希）
 * @returns {Promise<boolean>} 是否需要翻译
 */
async function needsTranslation(
  sourcePath: string,
  targetPath: string,
  state: TranslationState,
): Promise<boolean> {
  try {
    // 检查目标文件是否存在
    await fs.access(targetPath)
    // 内容哈希比对：只有源内容真正变化才重新翻译
    const source = await fs.readFile(sourcePath, 'utf-8')
    return state[toStateKey(sourcePath)] !== contentHash(source)
  } catch {
    return true
  }
}

interface TranslateFileOptions {
  force?: boolean
  skipExisting?: boolean
  outputSuffix?: string
}

/**
 * md 翻译后，用 zh→en 标签映射查表替换 frontmatter 的 tags（一对一、数量守恒，
 * en 界面显示英文标签）。缺失映射的标签保持中文并提示（tagMerger 会增量补齐映射）。
 * 支持单行数组（tags: ["a", "b"]）与多行数组（tags:\n  - a\n  - b）。
 */
let translationCache: Record<string, string> | null = null

function loadTranslationMapping(): Record<string, string> {
  if (translationCache) return translationCache
  translationCache = {}
  try {
    const mappingPath = path.join(cacheDir, 'tag-mapping.json')
    const raw = fsSync.readFileSync(mappingPath, 'utf-8')
    const parsed = JSON.parse(raw) as { translation?: Record<string, string> }
    if (parsed?.translation && typeof parsed.translation === 'object') {
      translationCache = parsed.translation
    }
  } catch (e) {
    console.warn(
      `[Warn] 标签映射加载失败（en 标签将保持中文）: ${e instanceof Error ? e.message : e}`,
    )
  }
  return translationCache
}

function translateFrontmatterTags(source: string, translated: string): string {
  const srcMatch = source.match(/^---\n([\s\S]*?)\n---/)
  if (!srcMatch) return translated
  const srcTags = srcMatch[1].match(/^tags:\s*(\[[\s\S]*?\])/m)
  if (!srcTags) return translated

  const translation = loadTranslationMapping()
  const zhTags = JSON.parse(srcTags[1]) as string[]
  const enTags = zhTags.map((tag) => {
    const en = translation[tag]
    if (!en) console.warn(`[Warn] 标签「${tag}」缺少 zh→en 映射（保持中文，请运行 tagMerger 补齐）`)
    return en || tag
  })
  const enTagsJson = JSON.stringify(enTags)
  // 只替换 LLM 输出的 frontmatter 区域，避免 tags 正则吞掉正文内容
  const fmMatch = translated.match(/^---\r?\n([\s\S]*?)\r?\n---/)
  if (!fmMatch) return translated
  const fm = fmMatch[1].replace(/^tags:\s*\[[\s\S]*?\]/m, `tags: ${enTagsJson}`)
  return `---\n${fm}\n---${translated.slice(fmMatch[0].length)}`
}

/**
 * 翻译单个文件
 * @param {string} inputFilePath - 输入文件路径
 * @param {Object} options - 选项
 * @returns {Promise<void>}
 */
async function translateFile(
  inputFilePath: string,
  options: TranslateFileOptions = {},
  state?: TranslationState,
): Promise<string | null> {
  const {
    force = false, // 强制重新翻译
    skipExisting = true, // 跳过已存在的翻译文件
    outputSuffix = '-en', // 输出文件后缀
  } = options

  // 独立调用时加载状态；目录批量调用由 translateDirectory 统一管理
  const dir = path.dirname(inputFilePath)
  const localState = state ?? loadState(dir)
  const saveLocal = !state

  try {
    // 检查文件扩展名
    const ext = path.extname(inputFilePath).toLowerCase()
    if (!['.md', '.yaml', '.yml'].includes(ext)) {
      console.log(`[WARN] 跳过不支持的文件类型: ${inputFilePath}`)
      return null
    }

    // 生成输出路径
    const basename = path.basename(inputFilePath, ext)
    const outputPath = path.join(path.dirname(inputFilePath), `${basename}${outputSuffix}${ext}`)

    // 检查是否需要翻译（内容哈希增量判定）
    if (skipExisting && !force) {
      const shouldTranslate = await needsTranslation(inputFilePath, outputPath, localState)
      if (!shouldTranslate) {
        console.log(`[INFO] 跳过已翻译文件: ${inputFilePath}`)
        return null // 返回null而不是undefined
      }
    }

    // 读取文件内容
    const content = await fs.readFile(inputFilePath, 'utf-8')

    // 翻译
    console.log(`[INFO] 正在翻译文件: ${inputFilePath}`)
    const fileType = ext === '.yaml' || ext === '.yml' ? 'yaml' : 'md'
    let translated = await translateText(content, fileType)

    // md 文章：frontmatter 的 tags 用 zh→en 映射查表翻译（一对一、数量守恒，
    // en 界面显示英文标签；缺失映射的标签保持中文并告警）
    if (fileType === 'md') {
      translated = translateFrontmatterTags(content, translated)
    }

    // 写入翻译后文件
    await fs.writeFile(outputPath, translated, 'utf-8')
    console.log(`[INFO] 翻译完成: ${outputPath}`)
    // 记录源内容哈希，下次跳过
    localState[toStateKey(inputFilePath)] = contentHash(content)
    if (saveLocal) saveState(dir, localState)

    return outputPath
  } catch (error) {
    console.error(`[ERROR] 处理文件时出错 (${inputFilePath}):`, (error as Error).message)
    throw error
  }
}

interface TranslateDirectoryOptions extends TranslateFileOptions {
  recursive?: boolean
  filePatterns?: string[]
  excludePatterns?: string[]
  concurrency?: number
}

/**
 * 批量翻译目录中的文件
 * @param {string} directoryPath - 目录路径
 * @param {Object} options - 选项
 * @returns {Promise<Array>} 翻译的文件列表
 */
async function translateDirectory(
  directoryPath: string,
  options: TranslateDirectoryOptions = {},
): Promise<string[]> {
  const {
    recursive = true, // 是否递归搜索子目录
    filePatterns = ['*.md', '*.yaml', '*.yml'], // 文件模式
    excludePatterns = ['*-en.*'], // 排除模式
    ...translateOptions
  } = options

  try {
    const files: string[] = []

    // 递归搜索文件
    async function searchFiles(dir: string) {
      const entries = await fs.readdir(dir, { withFileTypes: true })

      for (const entry of entries) {
        const fullPath = path.join(dir, entry.name)

        if (entry.isDirectory() && recursive) {
          await searchFiles(fullPath)
        } else if (entry.isFile()) {
          // 检查文件是否匹配模式且不被排除
          const matchesPattern = filePatterns.some((pattern) => {
            const regex = new RegExp(pattern.replace('*', '.*').replace('.', '\\.'))
            return regex.test(entry.name)
          })

          const isExcluded = excludePatterns.some((pattern) => {
            const regex = new RegExp(pattern.replace('*', '.*').replace('.', '\\.'))
            return regex.test(entry.name)
          })

          if (matchesPattern && !isExcluded) {
            files.push(fullPath)
          }
        }
      }
    }

    await searchFiles(directoryPath)

    console.log(`[INFO] 找到 ${files.length} 个需要翻译的文件`)

    // 统一加载翻译状态（内容哈希增量，杜绝 mtime 误判导致的重复翻译）
    const state = loadState(directoryPath)

    // 并发翻译：LLM 调用是 IO 密集，默认 4 路并发显著提速
    const concurrency = Math.max(1, Math.min(8, options.concurrency ?? 4))
    const results: string[] = []
    let cursor = 0
    const workers = Array.from({ length: concurrency }, async () => {
      while (cursor < files.length) {
        const file = files[cursor++]
        try {
          const result = await translateFile(file, translateOptions, state)
          if (result) results.push(result)
        } catch (error) {
          console.error(`[ERROR] 翻译失败: ${file}`, (error as Error).message)
        }
      }
    })
    await Promise.all(workers)

    // 保存翻译状态（下次仅翻译内容变化的源文件）
    saveState(directoryPath, state)

    console.log(`[INFO] 批量翻译完成，成功翻译 ${results.length} 个文件（并发 ${concurrency}）`)
    return results
  } catch (error) {
    console.error('[ERROR] 批量翻译时出错:', (error as Error).message)
    throw error
  }
}

/**
 * 翻译管理器 - 支持多种翻译策略
 */
class TranslationManager {
  options: TranslateDirectoryOptions

  constructor(options: TranslateDirectoryOptions = {}) {
    this.options = options
  }

  /**
   * 翻译单个文件或目录
   * @param {string} targetPath - 目标路径
   * @param {Object} options - 选项
   */
  async translate(
    targetPath: string,
    options: TranslateDirectoryOptions = {},
  ): Promise<string[] | string | null> {
    const mergedOptions = { ...this.options, ...options }

    try {
      const stats = await fs.stat(targetPath)

      if (stats.isDirectory()) {
        return await translateDirectory(targetPath, mergedOptions)
      } else {
        return await translateFile(targetPath, mergedOptions)
      }
    } catch (error) {
      console.error('[ERROR] 路径检查失败:', (error as Error).message)
      throw error
    }
  }

  /**
   * 强制重新翻译所有文件
   * @param {string} targetPath - 目标路径
   */
  async forceTranslate(targetPath: string) {
    return this.translate(targetPath, { force: true, skipExisting: false })
  }

  /**
   * 仅翻译新文件（跳过所有已存在的翻译）
   * @param {string} targetPath - 目标路径
   */
  async translateNewOnly(targetPath: string) {
    return this.translate(targetPath, { skipExisting: true })
  }
}

// 导出函数和类
export { translateText, translateFile, translateDirectory, TranslationManager, needsTranslation }

// 默认导出管理器实例
export default new TranslationManager()
