/**
 * Shared utilities for content generators.
 * Behavior must stay byte-identical to the original per-generator copies;
 * verify by diffing generated JSON before/after refactor.
 */

import fs from 'fs'
import path from 'path'
import yaml from 'js-yaml'

export function ensureDirectoryExistence(filePath: string) {
  const dirname = path.dirname(filePath)
  if (!fs.existsSync(dirname)) {
    fs.mkdirSync(dirname, { recursive: true })
  }
}

export function walk(
  dir: string,
  predicate: (filePath: string) => boolean = () => true,
  acc: string[] = [],
): string[] {
  if (!fs.existsSync(dir)) return acc
  const items = fs.readdirSync(dir, { withFileTypes: true })
  for (const it of items) {
    const full = path.join(dir, it.name)
    if (it.isDirectory()) walk(full, predicate, acc)
    else if (predicate(full)) acc.push(full)
  }
  return acc
}

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

// Crude markdown -> plain text for counting words
export function markdownToPlain(text: string): string {
  // remove code fences first
  let t = text.replace(/```[\s\S]*?```/g, ' ')
  // remove inline code
  t = t.replace(/`[^`]*`/g, ' ')
  // remove images ![alt](url)
  t = t.replace(/!\[[^\]]*]\([^)]+\)/g, ' ')
  // remove links [text](url)
  t = t.replace(/\[[^\]]*]\([^)]+\)/g, ' ')
  // remove html tags
  t = t.replace(/<[^>]+>/g, ' ')
  // remove headings/bold/italic/quotes/lists/tables markers
  t = t.replace(/^#{1,6}\s+/gm, ' ')
  t = t.replace(/[*_~`>#|-]{1,}/g, ' ')
  // remove reference links [id]: url
  t = t.replace(/^\s*\[[^\]]+]:\s+\S+.*$/gm, ' ')
  // collapse whitespace
  t = t.replace(/\s+/g, ' ').trim()
  return t
}

export function countWordsSmart(text: string): number {
  const cjk = (text.match(/[\u4E00-\u9FFF\u3400-\u4DBF]/g) || []).length
  const words = (text.match(/[A-Za-z0-9]+(?:'[A-Za-z0-9]+)?/g) || []).length
  return cjk + words
}

export function toPosixRelativeNoExt(fullPath: string, baseDir: string): string {
  const rel = path.relative(baseDir, fullPath)
  const noExt = rel.replace(/\.[^/.]+$/, '')
  return noExt.split(path.sep).join('/')
}

export function extractTitleFromH1(raw: string): string {
  const lines = raw.split(/\r?\n/)
  for (const line of lines) {
    const m = line.match(/^#\s+(.*)$/)
    if (m) return m[1].trim()
  }
  return ''
}
