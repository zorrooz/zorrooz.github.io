/**
 * TagMerger — AI 辅助标签合并工具（核心逻辑）
 *
 * 收集博客全部标签（notes md frontmatter + categories.yaml 的 projects/topics tags），
 * 调用 LLM 生成「旧标签 → 规范标签」的合并映射，支持 dry-run 预览与 apply 写回。
 */

import fs from 'fs'
import path from 'path'
import yaml from 'js-yaml'

import { contentSrcDir, cacheDir, enSrcDir, localeSuffix, EN_SUFFIX } from '../dataConfig.ts'
import { walk, parseFrontMatterAndBody, normalizeTags } from '../generators/core/index.ts'
import { completeChat, loadLlmConfig, type LlmConfig } from '../lib/llm.ts'
import { rewriteFrontmatterTags } from '../lib/frontmatter.ts'
import { countTags } from '../lib/tags.ts'
import { loadTagMapping, type TagMapping } from '../lib/tagMapping.ts'

export { type TagMapping }

export const defaultMappingPath = path.join(cacheDir, 'tag-mapping.json')

/** 收集指定 locale 下全部标签（md frontmatter + categories yaml），返回 标签 → 出现次数 */
export function collectTags(locale: 'zh-CN' | 'en-US'): Map<string, number> {
  const map = new Map<string, number>()

  const addTags = (tags: unknown) => {
    for (const [t, count] of countTags(normalizeTags(tags))) {
      map.set(t, (map.get(t) || 0) + count)
    }
  }

  // 目录即契约：zh 扫手写源，en 扫机器翻译层（各自目录内全部 .md）
  const srcDir = locale === 'en-US' ? enSrcDir : contentSrcDir
  const mdFiles = walk(srcDir, (p) => /\.md$/i.test(p))
  for (const mdPath of mdFiles) {
    const raw = fs.readFileSync(mdPath, 'utf-8')
    addTags(parseFrontMatterAndBody(raw).frontmatter.tags)
  }

  const yamlPath = path.join(srcDir, `categories${localeSuffix(locale)}.yaml`)
  if (fs.existsSync(yamlPath)) {
    try {
      const obj = yaml.load(fs.readFileSync(yamlPath, 'utf-8')) as Record<string, unknown>
      for (const key of ['projects', 'topics']) {
        const arr = Array.isArray(obj?.[key]) ? (obj[key] as Record<string, unknown>[]) : []
        for (const item of arr) addTags(item?.tags)
      }
    } catch (e) {
      console.warn(`Warn: failed to parse ${yamlPath}:`, e instanceof Error ? e.message : e)
    }
  }

  return map
}

/** 将标签统计转为 LLM 输入文本（按次数降序） */
export function formatTagStats(tags: Map<string, number>): string {
  const sorted = Array.from(tags.entries())
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], 'zh-Hans-CN'))
    .map(([tag, count]) => `${tag} (${count})`)
  return sorted.join('\n')
}

/** 解析 LLM 输出 JSON：剥代码围栏 → 取首尾大括号 → JSON.parse → 校验对象 → 过滤非字符串，取指定键 */
export function extractJsonObject(text: string, key: string): Record<string, string> {
  let t = text.trim()
  const fence = t.match(/```(?:json)?\s*([\s\S]*?)```/i)
  if (fence) t = fence[1].trim()
  const start = t.indexOf('{')
  const end = t.lastIndexOf('}')
  if (start === -1 || end <= start) throw new Error('LLM 输出中未找到 JSON')
  const parsed = JSON.parse(t.slice(start, end + 1)) as Record<string, unknown>
  const obj = parsed[key]
  if (!obj || typeof obj !== 'object' || Array.isArray(obj)) {
    throw new Error(`LLM 输出缺少 ${key} 对象`)
  }
  const result: Record<string, string> = {}
  for (const [k, v] of Object.entries(obj as Record<string, unknown>)) {
    if (typeof v === 'string' && v.trim()) result[k] = v.trim()
  }
  return result
}

/** 解析 LLM 输出的标签合并映射 JSON（{"mapping": {...}}） */
export function extractJson(text: string): Record<string, string> {
  return extractJsonObject(text, 'mapping')
}

/** 解析 LLM 输出的 zh→en 翻译映射 JSON（{"translation": {...}}） */
export function extractTranslationJson(text: string): Record<string, string> {
  return extractJsonObject(text, 'translation')
}

/** 找出 translation 映射中缺失的 zh 标签（增量翻译只用这部分，省 token） */
export function missingTranslationTags(
  zhTags: Map<string, number>,
  translation: Record<string, string>,
): Map<string, number> {
  const missing = new Map<string, number>()
  for (const [tag, count] of zhTags) {
    if (!translation[tag]) missing.set(tag, count)
  }
  return missing
}

/** LLM 生成 zh→en 标签翻译映射（一对一、数量守恒；仅传缺失标签） */
export async function generateTranslationWithLlm(
  zhTags: Map<string, number>,
  api: LlmConfig,
): Promise<Record<string, string>> {
  const systemPrompt = `你是专业翻译。下面是一个中文技术博客的标签列表（括号内为出现次数）。
请把每个标签翻译为英文，要求：
1. 输出 JSON：{"translation": {"中文标签": "英文标签"}}
2. 一对一映射，标签数量与输入完全一致
3. 专有名词/软件名保持原样（如 Python、BWA、HISAT2、ggplot2、KaTeX、RELION、PyMOL、ChimeraX）
4. 英文每个单词首字母大写、其余小写
5. 只输出 JSON，不要任何解释或额外内容`
  const result = await completeChat(systemPrompt, formatTagStats(zhTags), {
    temperature: 0.2,
    thinking: api.thinking,
  })
  if (!result.text) throw new Error('LLM 返回空内容')
  return extractTranslationJson(result.text)
}

/**
 * 确保 zh→en 标签翻译映射存在（构建时自动调用）：
 * 缓存命中（所有 zh 标签已有映射）→ 0 token；仅缺失标签增量生成。
 * 映射写入 tag-mapping.json 的 translation 段。
 */
export async function ensureTagTranslation(
  mappingPath: string = defaultMappingPath,
): Promise<Record<string, string>> {
  const zhTags = collectTags('zh-CN')
  const existing = loadMapping(mappingPath)
  const translation: Record<string, string> = existing?.translation || {}

  const missing = missingTranslationTags(zhTags, translation)
  if (missing.size === 0) {
    return translation
  }

  const llm = await loadLlmConfig()
  if (!llm) {
    console.warn('[Warn] tagMerger: 未找到 LLM 配置，无法生成 zh→en 标签映射（en 标签将保持中文）')
    return translation
  }

  console.log(`[INFO] tagMerger: 生成 ${missing.size} 个缺失标签的 zh→en 映射...`)
  const generated = await generateTranslationWithLlm(missing, llm)
  const merged = { ...translation, ...generated }
  const out: TagMapping = {
    locale: 'zh-CN',
    generatedAt: new Date().toISOString(),
    mapping: existing?.mapping || {},
    translation: merged,
  }
  fs.writeFileSync(mappingPath, JSON.stringify(out, null, 2), 'utf-8')
  console.log(`[INFO] tagMerger: 映射已更新 ${Object.keys(generated).length} 项 → ${mappingPath}`)
  return merged
}

/**
 * 自动解决中英标签不一致（构建时调用）：
 * 以 zh 文件 frontmatter/yaml 的 tags 为唯一基准，经 zh→en 映射翻译后重写 cache/en 的 -en 文件。
 * 覆盖：notes 的 -en.md 与 categories-en.yaml（projects/topics 的 tags）。
 * 中英文件按相对结构配对：cache/en 镜像 content-src（-en.md ↔ .md）。
 * 缺映射的标签保持中文并告警（ensureTagTranslation 会增量补齐映射）。
 */
export function fixTagConsistency(): { fixed: number; missing: number } {
  const mapping = loadMapping(defaultMappingPath)
  const translation = mapping?.translation || {}
  let fixed = 0
  let missing = 0
  const missingTags = new Set<string>()

  const translateList = (zhTags: string[]): string[] =>
    zhTags.map((t) => {
      const en = translation[t]
      if (!en) {
        missing++
        missingTags.add(t)
        return t
      }
      return en
    })

  // 1. cache/en 的 -en.md：以镜像位置对应的 content-src zh.md 为基准
  for (const enPath of walk(enSrcDir, (p) => /\.md$/i.test(p))) {
    const rel = path.relative(enSrcDir, enPath)
    const zhPath = path.join(contentSrcDir, rel.replace(new RegExp(`${EN_SUFFIX}\\.md$`), '.md'))
    if (!fs.existsSync(zhPath)) continue
    const zh = parseFrontMatterAndBody(fs.readFileSync(zhPath, 'utf-8'))
    if (!Array.isArray(zh.frontmatter.tags)) continue
    const enTags = translateList(zh.frontmatter.tags as string[])
    const raw = fs.readFileSync(enPath, 'utf-8')
    // 只替换 frontmatter 区域的 tags 行/块（见 lib/frontmatter.ts），避免正则吞掉正文内容
    const updated = rewriteFrontmatterTags(raw, () => enTags)
    if (updated !== raw) {
      fs.writeFileSync(enPath, updated, 'utf-8')
      fixed++
      console.log(`[FIX] 重写 tags: ${rel} -> ${enTags.join(', ')}`)
    }
  }

  // 2. categories-en.yaml：以 categories.yaml 为基准（projects/topics 的 tags）
  //    按行序配对（zh/en yaml 结构由翻译器同步维护，顺序一致）——不能用「缩进+tags:」做
  //    Map key，同一缩进的多个 tags 行会互相覆盖。
  const zhYamlPath = path.join(contentSrcDir, 'categories.yaml')
  const enYamlPath = path.join(enSrcDir, `categories${EN_SUFFIX}.yaml`)
  if (fs.existsSync(zhYamlPath) && fs.existsSync(enYamlPath)) {
    const zhYaml = fs.readFileSync(zhYamlPath, 'utf-8')
    const enYaml = fs.readFileSync(enYamlPath, 'utf-8')
    const zhTagValues: string[] = []
    for (const line of zhYaml.split('\n')) {
      const m = line.match(/^(\s*tags:\s*)(\[.*\])$/)
      if (m) {
        try {
          zhTagValues.push(JSON.stringify(translateList(JSON.parse(m[2]) as string[])))
        } catch {
          zhTagValues.push(m[2])
        }
      }
    }
    if (zhTagValues.length) {
      const lines = enYaml.split('\n')
      let k = 0
      let changed = false
      for (let i = 0; i < lines.length; i++) {
        const m = lines[i].match(/^(\s*tags:\s*)(\[.*\])$/)
        if (m && k < zhTagValues.length) {
          const fixedLine = `${m[1]}${zhTagValues[k]}`
          if (lines[i] !== fixedLine) {
            lines[i] = fixedLine
            changed = true
          }
          k++
        }
      }
      if (changed) {
        fs.writeFileSync(enYamlPath, lines.join('\n'), 'utf-8')
        fixed++
        console.log('[FIX] 重写 tags: categories-en.yaml')
      }
    }
  }

  if (missingTags.size) {
    console.warn(
      `[Warn] ${missingTags.size} 个标签缺 zh→en 映射（保持中文，构建会自动补齐映射）: ${[
        ...missingTags,
      ].join(', ')}`,
    )
  }
  console.log(`[tag-consistency] 修复 ${fixed} 个文件，缺映射 ${missing} 个`)
  return { fixed, missing }
}

export async function generateMappingWithLlm(
  tags: Map<string, number>,
  locale: 'zh-CN' | 'en-US',
  api: LlmConfig,
): Promise<Record<string, string>> {
  const isZh = locale === 'zh-CN'
  const systemPrompt = isZh
    ? `你是标签规范化专家。下面是一个中文技术博客的全部标签（括号内为使用次数）。
请输出标签合并映射 JSON，规则：
1. 合并语义相同或相近的标签（如 "生物信息" 与 "生物信息学"、"R" 与 "R语言"）
2. 专有名词（软件名、工具名如 BWA、HISAT2、Pandas、Biopython、MaxQuant、ggplot2）必须保留原样
3. 规范标签优先选择更通用、更常用的那个
4. 不需要合并的标签也列出，映射到自身
5. 只输出 JSON：{"mapping": {"旧标签": "新标签"}}，不要任何解释`
    : `You are a tag normalization expert. Below is the full tag list of an English technical blog (usage count in parentheses).
Output a tag merge mapping JSON following these rules:
1. Merge semantically identical or similar tags (e.g. "RNA-seq" and "RNA Seq", "R" and "R language")
2. Proper nouns (software/tool names such as BWA, HISAT2, Pandas, Biopython, MaxQuant, ggplot2) must be kept unchanged
3. Prefer the more general / more common tag as the canonical one
4. List every tag, mapping unchanged tags to themselves
5. Output only JSON: {"mapping": {"old tag": "new tag"}}, no explanations`

  const result = await completeChat(systemPrompt, formatTagStats(tags), {
    temperature: 0.2,
    thinking: api.thinking,
  })
  if (!result.text) throw new Error('LLM 返回空内容')
  return extractJson(result.text)
}

/** 读取已存在的映射文件（完整逻辑见 lib/tagMapping.ts） */
export function loadMapping(mappingPath: string): TagMapping | null {
  return loadTagMapping(mappingPath)
}

/** 将映射写回 md frontmatter 的 tags 字段（仅替换 tags 行/块，保留其余格式与行尾） */
export function applyMappingToMarkdown(
  mapping: Record<string, string>,
  locale: 'zh-CN' | 'en-US' = 'zh-CN',
): number {
  let changed = 0
  const srcDir = locale === 'en-US' ? enSrcDir : contentSrcDir
  const mdFiles = walk(srcDir, (p) => /\.md$/i.test(p))
  for (const mdPath of mdFiles) {
    const raw = fs.readFileSync(mdPath, 'utf-8')
    const { frontmatter } = parseFrontMatterAndBody(raw)
    if (!frontmatter.tags) continue

    const oldTags = normalizeTags(frontmatter.tags)
    if (!oldTags.length) continue
    const newTags = Array.from(new Set(oldTags.map((t) => mapping[t] || t)))
    if (newTags.join('\u0000') === oldTags.join('\u0000')) continue
    if (newTags.some((t) => /[,"[\]\n\r]/.test(t))) {
      console.warn(`Warn: 标签含特殊字符（逗号/引号/方括号/换行），跳过文件：${mdPath}`)
      continue
    }

    const updated = rewriteFrontmatterTags(raw, () => newTags)
    if (updated === raw) continue
    fs.writeFileSync(mdPath, updated, 'utf-8')
    changed++
  }
  return changed
}

/** 将映射写回 categories.yaml 的 tags 字段（行内数组替换，保留文件其余格式） */
export function applyMappingToCategoriesYaml(
  mapping: Record<string, string>,
  locale: 'zh-CN' | 'en-US' = 'zh-CN',
): number {
  const yamlPath = path.join(
    locale === 'en-US' ? enSrcDir : contentSrcDir,
    `categories${localeSuffix(locale)}.yaml`,
  )
  if (!fs.existsSync(yamlPath)) return 0

  const raw = fs.readFileSync(yamlPath, 'utf-8')
  const eol = raw.includes('\r\n') ? '\r\n' : '\n'
  const lines = raw.split(/\r?\n/)
  let changed = false
  const out = lines.map((line) => {
    const m = line.match(/^(\s*tags:\s*\[)([^\]]*)(\])$/)
    if (!m) return line
    const tags = m[2]
      .split(',')
      .map((s) => s.trim().replace(/^"|"$/g, ''))
      .filter(Boolean)
    const newTags = Array.from(new Set(tags.map((t) => mapping[t] || t)))
    if (newTags.join(',') === tags.join(',')) return line
    if (newTags.some((t) => /[,"[\]\n\r]/.test(t))) {
      console.warn(`Warn: 标签含特殊字符（逗号/引号/方括号/换行），跳过该 tags 行：${line}`)
      return line
    }
    changed = true
    return `${m[1]}${newTags.map((t) => JSON.stringify(t)).join(', ')}${m[3]}`
  })
  if (changed) {
    fs.writeFileSync(yamlPath, out.join(eol), 'utf-8')
    return 1
  }
  return 0
}

export { loadLlmConfig }
