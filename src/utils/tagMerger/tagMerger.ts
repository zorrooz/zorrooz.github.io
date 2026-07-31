/**
 * TagMerger — AI 辅助标签合并工具（核心逻辑）
 *
 * 收集博客全部标签（notes md frontmatter + categories.yaml 的 projects/topics tags），
 * 调用 LLM 生成「旧标签 → 规范标签」的合并映射，支持 dry-run 预览与 apply 写回。
 */

import fs from 'fs'
import path from 'path'
import { fileURLToPath, pathToFileURL } from 'url'
import yaml from 'js-yaml'
import OpenAI from 'openai'

import { walk, parseFrontMatterAndBody, normalizeTags } from '../generators/core/index.ts'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

export const contentSrcDir = path.join(__dirname, '../../content-src')
export const defaultMappingPath = path.join(contentSrcDir, 'tag-mapping.json')

export interface TagMapping {
  locale: string
  generatedAt: string
  mapping: Record<string, string>
}

export interface TagStat {
  tag: string
  count: number
}

/** 收集指定 locale 下全部标签（md frontmatter + categories yaml），返回 标签 → 出现次数 */
export function collectTags(locale: 'zh-CN' | 'en-US'): Map<string, number> {
  const map = new Map<string, number>()
  const suffix = locale === 'en-US' ? '-en' : ''

  const addTags = (tags: unknown) => {
    for (const t of normalizeTags(tags)) map.set(t, (map.get(t) || 0) + 1)
  }

  const mdFiles = walk(contentSrcDir, (p) => {
    if (!/\.md$/i.test(p)) return false
    return locale === 'en-US' ? p.endsWith('-en.md') : !p.endsWith('-en.md')
  })
  for (const mdPath of mdFiles) {
    const raw = fs.readFileSync(mdPath, 'utf-8')
    addTags(parseFrontMatterAndBody(raw).frontmatter.tags)
  }

  const yamlPath = path.join(contentSrcDir, `categories${suffix}.yaml`)
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

export function extractJson(text: string): Record<string, string> {
  let t = text.trim()
  const fence = t.match(/```(?:json)?\s*([\s\S]*?)```/i)
  if (fence) t = fence[1].trim()
  const start = t.indexOf('{')
  const end = t.lastIndexOf('}')
  if (start === -1 || end <= start) throw new Error('LLM 输出中未找到 JSON')
  const parsed = JSON.parse(t.slice(start, end + 1)) as Record<string, unknown>
  const mapping = parsed.mapping
  if (!mapping || typeof mapping !== 'object' || Array.isArray(mapping)) {
    throw new Error('LLM 输出缺少 mapping 对象')
  }
  const result: Record<string, string> = {}
  for (const [oldTag, newTag] of Object.entries(mapping as Record<string, unknown>)) {
    if (typeof newTag === 'string' && newTag.trim()) result[oldTag] = newTag.trim()
  }
  return result
}

export async function generateMappingWithLlm(
  tags: Map<string, number>,
  locale: 'zh-CN' | 'en-US',
  api: { url: string; apikey: string; model: string },
): Promise<Record<string, string>> {
  const openai = new OpenAI({ baseURL: api.url, apiKey: api.apikey })

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

  const completion = await openai.chat.completions.create({
    messages: [
      { role: 'system', content: systemPrompt },
      { role: 'user', content: formatTagStats(tags) },
    ],
    model: api.model,
    temperature: 0.2,
  })

  const content = completion.choices[0]?.message?.content
  if (!content) throw new Error('LLM 返回空内容')
  return extractJson(content)
}

async function loadLlmConfig(): Promise<{ url: string; apikey: string; model: string } | null> {
  const cfgPath = path.join(__dirname, '../../config/llmConfig.ts')
  if (!fs.existsSync(cfgPath)) return null
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

/** 读取已存在的映射文件 */
export function loadMapping(mappingPath: string): TagMapping | null {
  try {
    if (!fs.existsSync(mappingPath)) return null
    const parsed = JSON.parse(fs.readFileSync(mappingPath, 'utf-8')) as TagMapping
    if (!parsed?.mapping || typeof parsed.mapping !== 'object') return null
    return parsed
  } catch {
    return null
  }
}

/** 将映射写回 md frontmatter 的 tags 字段（仅替换 tags 行/块，保留其余格式与行尾） */
export function applyMappingToMarkdown(
  mapping: Record<string, string>,
  locale: 'zh-CN' | 'en-US' = 'zh-CN',
): number {
  let changed = 0
  const mdFiles = walk(contentSrcDir, (p) => {
    if (!/\.md$/i.test(p)) return false
    return locale === 'en-US' ? p.endsWith('-en.md') : !p.endsWith('-en.md')
  })
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

    const eol = raw.includes('\r\n') ? '\r\n' : '\n'
    const lines = raw.split(/\r?\n/)
    let tagsStart = -1
    let tagsEnd = -1
    for (let i = 0; i < lines.length; i++) {
      if (/^tags:/i.test(lines[i])) {
        tagsStart = i
        if (/^tags:\s*\[/i.test(lines[i])) {
          tagsEnd = i
        } else {
          let j = i + 1
          while (j < lines.length && /^\s+-\s+/.test(lines[j])) j++
          tagsEnd = j - 1
        }
        break
      }
    }
    if (tagsStart === -1) continue

    const rendered = newTags.map((t) => `  - ${JSON.stringify(t)}`).join(eol)
    const replacement =
      tagsEnd === tagsStart
        ? `tags: [${newTags.map((t) => JSON.stringify(t)).join(', ')}]`
        : `tags:${eol}${rendered}`
    lines.splice(tagsStart, tagsEnd - tagsStart + 1, replacement)
    fs.writeFileSync(mdPath, lines.join(eol), 'utf-8')
    changed++
  }
  return changed
}

/** 将映射写回 categories.yaml 的 tags 字段（行内数组替换，保留文件其余格式） */
export function applyMappingToCategoriesYaml(
  mapping: Record<string, string>,
  locale: 'zh-CN' | 'en-US' = 'zh-CN',
): number {
  const yamlPath = path.join(contentSrcDir, locale === 'en-US' ? 'categories-en.yaml' : 'categories.yaml')
  if (!fs.existsSync(yamlPath)) return 0

  const raw = fs.readFileSync(yamlPath, 'utf-8')
  const eol = raw.includes('\r\n') ? '\r\n' : '\n'
  const lines = raw.split(/\r?\n/)
  let changed = false
  const out = lines.map((line) => {
    const m = line.match(/^(\s*tags:\s*\[)([^\]]*)(\])$/)
    if (!m) return line
    const tags = m[2].split(',').map((s) => s.trim().replace(/^"|"$/g, '')).filter(Boolean)
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
