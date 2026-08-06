import OpenAI from 'openai'
import fs from 'fs/promises'
import fsSync from 'fs'
import path from 'path'
import crypto from 'crypto'
import { pathToFileURL } from 'url'
import { cacheDir, contentSrcDir, enSrcDir } from '../dataConfig.ts'

let clientCache: OpenAI | null = null

/**
 * 懒加载 DeepSeek 配置（llmConfig.ts 被 gitignore，CI 无此文件）。
 * 与 tagMerger 相同：模块可能被 vite-ssg 打包，静态 import 会因缺文件失败。
 */
async function loadLlmConfig(): Promise<{ url: string; apikey: string; model: string } | null> {
  const cfgPath = path.join(import.meta.dirname, '../llmConfig.ts')
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
  if (!cfg)
    throw new Error('未找到 LLM 配置：请在 scripts/llmConfig.ts 配置 { url, apikey, model }')
  clientCache = new OpenAI({ baseURL: cfg.url, apiKey: cfg.apikey })
  return clientCache
}

/** 翻译一段文本（md 或 yaml），返回英文结果 */
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

/** 内容哈希（增量判定依据）：比 mtime 可靠，复制/同步/换分支不误触翻译 */
function contentHash(content: string): string {
  return crypto.createHash('md5').update(content, 'utf-8').digest('hex')
}

/** 状态键：相对 content-src 的路径，使 .translate-state.json 可跨机器复用 */
function toStateKey(sourcePath: string): string {
  return path.relative(contentSrcDir, sourcePath).split(path.sep).join('/')
}

interface TranslationState {
  [sourcePath: string]: string
}

/** 翻译状态统一存 cache/（数据分支第二层），键为相对 content-src 的路径 */
function stateFilePath(): string {
  return path.join(cacheDir, '.translate-state.json')
}

function loadState(): TranslationState {
  try {
    return JSON.parse(fsSync.readFileSync(stateFilePath(), 'utf-8')) as TranslationState
  } catch {
    return {}
  }
}

function saveState(state: TranslationState): void {
  // 空状态不写，避免覆盖已建立的翻译状态（否则下次全量重翻）
  if (!state || Object.keys(state).length === 0) return
  try {
    fsSync.writeFileSync(stateFilePath(), JSON.stringify(state, null, 2), 'utf-8')
  } catch (e) {
    console.warn(
      `[Warn] 翻译状态写入失败（不影响本次翻译，但下次会重复翻译）: ${
        e instanceof Error ? e.message : e
      }`,
    )
  }
}

/** 目标文件缺失或源内容哈希变化时需要重新翻译 */
async function needsTranslation(
  sourcePath: string,
  targetPath: string,
  state: TranslationState,
): Promise<boolean> {
  try {
    await fs.access(targetPath)
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

/** md 翻译后按 zh→en 标签映射查表替换 frontmatter tags（数量守恒，缺失的保持中文并告警） */
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
 * 翻译单个文件。目标不存在或内容变化时才翻译（见 needsTranslation）；
 * 输出镜像 content-src 的相对结构写入 cache/en（第二层机器区，-en 为内容身份后缀）。
 */
async function translateFile(
  inputFilePath: string,
  options: TranslateFileOptions = {},
  state?: TranslationState,
): Promise<string | null> {
  const { force = false, skipExisting = true, outputSuffix = '-en' } = options

  const localState = state ?? loadState()
  const saveLocal = !state

  try {
    const ext = path.extname(inputFilePath).toLowerCase()
    if (!['.md', '.yaml', '.yml'].includes(ext)) {
      console.log(`[WARN] 跳过不支持的文件类型: ${inputFilePath}`)
      return null
    }

    const basename = path.basename(inputFilePath, ext)
    const rel = path.relative(contentSrcDir, inputFilePath)
    const outputPath =
      !rel.startsWith('..') && !path.isAbsolute(rel)
        ? path.join(enSrcDir, path.dirname(rel), `${basename}${outputSuffix}${ext}`)
        : path.join(path.dirname(inputFilePath), `${basename}${outputSuffix}${ext}`)

    if (skipExisting && !force) {
      const shouldTranslate = await needsTranslation(inputFilePath, outputPath, localState)
      if (!shouldTranslate) {
        console.log(`[INFO] 跳过已翻译文件: ${inputFilePath}`)
        return null
      }
    }

    const content = await fs.readFile(inputFilePath, 'utf-8')

    console.log(`[INFO] 正在翻译文件: ${inputFilePath}`)
    const fileType = ext === '.yaml' || ext === '.yml' ? 'yaml' : 'md'
    let translated = await translateText(content, fileType)

    if (fileType === 'md') {
      translated = translateFrontmatterTags(content, translated)
    }

    await fs.mkdir(path.dirname(outputPath), { recursive: true })
    await fs.writeFile(outputPath, translated, 'utf-8')
    console.log(`[INFO] 翻译完成: ${outputPath}`)
    localState[toStateKey(inputFilePath)] = contentHash(content)
    if (saveLocal) saveState(localState)

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

/** 批量翻译目录中的文件，默认 4 路并发调用 LLM，返回成功翻译的文件列表 */
async function translateDirectory(
  directoryPath: string,
  options: TranslateDirectoryOptions = {},
): Promise<string[]> {
  const {
    recursive = true,
    filePatterns = ['*.md', '*.yaml', '*.yml'],
    excludePatterns = ['*-en.*'],
    ...translateOptions
  } = options

  try {
    const files: string[] = []

    async function searchFiles(dir: string) {
      const entries = await fs.readdir(dir, { withFileTypes: true })

      for (const entry of entries) {
        const fullPath = path.join(dir, entry.name)

        if (entry.isDirectory() && recursive) {
          await searchFiles(fullPath)
        } else if (entry.isFile()) {
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

    const state = loadState()

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

    saveState(state)

    console.log(`[INFO] 批量翻译完成，成功翻译 ${results.length} 个文件（并发 ${concurrency}）`)
    return results
  } catch (error) {
    console.error('[ERROR] 批量翻译时出错:', (error as Error).message)
    throw error
  }
}

/** 翻译管理器：面向单个文件或目录的统一入口 */
class TranslationManager {
  options: TranslateDirectoryOptions

  constructor(options: TranslateDirectoryOptions = {}) {
    this.options = options
  }

  /** 翻译单个文件或目录 */
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

  /** 强制重新翻译所有文件 */
  async forceTranslate(targetPath: string) {
    return this.translate(targetPath, { force: true, skipExisting: false })
  }

  /** 仅翻译新文件（跳过所有已存在的翻译） */
  async translateNewOnly(targetPath: string) {
    return this.translate(targetPath, { skipExisting: true })
  }
}

export { translateText, translateFile, translateDirectory, TranslationManager, needsTranslation }

export default new TranslationManager()
