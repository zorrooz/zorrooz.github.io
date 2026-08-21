// scripts/lib/frontmatter.ts — frontmatter 解析与 tags 改写共享实现（translator 与 tagMerger 共用）
import yaml from 'js-yaml'

export function normalizeTags(tags: unknown): string[] {
  if (Array.isArray(tags)) {
    return tags.map((t) => (typeof t === 'string' ? t.trim() : '')).filter(Boolean)
  }
  if (typeof tags === 'string') {
    return tags
      .split(/[,，]/g)
      .map((s) => s.trim())
      .filter(Boolean)
  }
  return []
}

interface FrontMatterResult {
  frontmatter: Record<string, unknown>
  body: string
}

export function parseFrontMatterAndBody(raw: string): FrontMatterResult {
  if (raw.startsWith('---')) {
    const lines = raw.split(/\r?\n/)
    let endIdx = -1
    for (let i = 1; i < lines.length; i++) {
      if (lines[i].trim() === '---') {
        endIdx = i
        break
      }
    }
    if (endIdx !== -1) {
      const fmText = lines.slice(1, endIdx).join('\n')
      let fm: Record<string, unknown> = {}
      try {
        fm = (yaml.load(fmText) as Record<string, unknown>) || {}
      } catch (e) {
        console.warn(
          'Warn: failed to parse frontmatter YAML. Using empty object.',
          e instanceof Error ? e.message : e,
        )
      }
      const body = lines.slice(endIdx + 1).join('\n')
      return { frontmatter: fm || {}, body }
    }
  }
  return { frontmatter: {}, body: raw }
}

/** 解析 frontmatter（复用 parseFrontMatterAndBody 语义）；无 frontmatter 返回 null */
export function parseFrontmatter(raw: string): {
  data: Record<string, unknown>
  body: string
} | null {
  if (!raw.startsWith('---')) return null
  const { frontmatter, body } = parseFrontMatterAndBody(raw)
  if (body === raw) return null
  return { data: frontmatter, body }
}

/** 解析行内数组 `tags: [...]`（可能跨行）；非 JSON（裸值）时按逗号切分并剥引号 */
function parseInlineTags(span: string): string[] {
  const start = span.indexOf('[')
  if (start === -1) return []
  const end = span.lastIndexOf(']')
  const inner = (end > start ? span.slice(start + 1, end) : span.slice(start + 1)).trim()
  if (!inner) return []
  try {
    const parsed = JSON.parse(`[${inner}]`) as unknown
    if (Array.isArray(parsed)) {
      return parsed
        .filter((t): t is string => typeof t === 'string')
        .map((t) => t.trim())
        .filter(Boolean)
    }
  } catch {
    // 非 JSON（未加引号的裸值），走逗号切分
  }
  return inner
    .split(',')
    .map((s) => s.trim().replace(/^["']|["']$/g, ''))
    .filter(Boolean)
}

/** 解析块状 list：`tags:` 后跟若干 `  - item` 行 */
function parseBlockTags(lines: string[], start: number, end: number): string[] {
  const tags: string[] = []
  for (let i = start + 1; i <= end; i++) {
    const m = lines[i].match(/^\s+-\s+(.*)$/)
    if (!m) continue
    const t = m[1].trim().replace(/^["']|["']$/g, '')
    if (t) tags.push(t)
  }
  return tags
}

/**
 * 统一改写 frontmatter 的 tags 行/块：支持行内数组 `tags: [a, b]` 与块状 list
 * `tags:\n  - a` 两种语法，保留其余格式与行尾；无 frontmatter 或无 tags 时原样返回。
 */
export function rewriteFrontmatterTags(
  raw: string,
  mapFn: (tags: string[]) => string[],
): string {
  const lines = raw.split(/\r?\n/)
  if (lines[0]?.trim() !== '---') return raw
  let fmEnd = -1
  for (let i = 1; i < lines.length; i++) {
    if (lines[i].trim() === '---') {
      fmEnd = i
      break
    }
  }
  if (fmEnd === -1) return raw

  let tagsStart = -1
  let tagsEnd = -1
  let inline = false
  for (let i = 1; i < fmEnd; i++) {
    if (!/^tags:/.test(lines[i])) continue
    tagsStart = i
    if (/^tags:\s*\[/.test(lines[i])) {
      let j = i
      while (j < fmEnd && !lines[j].includes(']')) j++
      tagsEnd = j < fmEnd ? j : i
      inline = true
    } else {
      let j = i + 1
      while (j < fmEnd && /^\s+-\s+/.test(lines[j])) j++
      tagsEnd = j - 1
    }
    break
  }
  if (tagsStart === -1) return raw

  const eol = raw.includes('\r\n') ? '\r\n' : '\n'
  const oldTags = inline
    ? parseInlineTags(lines.slice(tagsStart, tagsEnd + 1).join('\n'))
    : parseBlockTags(lines, tagsStart, tagsEnd)
  const newTags = mapFn(oldTags)

  const replacement = inline
    ? `tags: [${newTags.map((t) => JSON.stringify(t)).join(', ')}]`
    : `tags:${eol}${newTags.map((t) => `  - ${JSON.stringify(t)}`).join(eol)}`
  const original = lines.slice(tagsStart, tagsEnd + 1).join(eol)
  if (replacement === original) return raw
  return [...lines.slice(0, tagsStart), replacement, ...lines.slice(tagsEnd + 1)].join(eol)
}
