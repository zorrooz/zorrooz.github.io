/**
 * Shared utilities for content generators.
 * Behavior must stay byte-identical to the original per-generator copies;
 * verify by diffing generated JSON before/after refactor.
 */

import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'
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

export function safeArray(x: unknown): unknown[] {
  return Array.isArray(x) ? x : []
}

/**
 * Parse a JSON file. Returns the parsed value (could be null when the file
 * literally contains `null`), or null when the file is missing / unreadable.
 * Callers distinguish "missing" via fs.existsSync; invalid JSON yields null.
 */
export function readJson(filePath: string): unknown {
  try {
    if (!fs.existsSync(filePath)) return null
    return JSON.parse(fs.readFileSync(filePath, 'utf-8').trim())
  } catch {
    return null
  }
}

/** Parse a YAML file; returns null when missing, or undefined for an empty file. */
export function readYaml(filePath: string): unknown {
  try {
    if (!fs.existsSync(filePath)) return null
    return yaml.load(fs.readFileSync(filePath, 'utf-8'))
  } catch {
    return null
  }
}

export function writeJsonFile(filePath: string, data: unknown, space = 2): void {
  ensureDirectoryExistence(filePath)
  fs.writeFileSync(filePath, JSON.stringify(data, null, space), 'utf-8')
}

/**
 * True when the module is executed directly (node foo.ts), tolerant of
 * relative/absolute argv forms. Used to gate CLI entry points so that
 * importing a generator module does not trigger its main().
 */
export function isDirectRun(importMeta: ImportMeta): boolean {
  return Boolean(
    process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(importMeta.url),
  )
}

/** Standard single-generator CLI bootstrap: banner + run + exit code on failure. */
export function runCliScript(name: string, run: () => void | Promise<void>): void {
  console.log(`Starting ${name} generation script...`)
  const done = () => console.log(`${name} generation complete.`)
  try {
    const result = run()
    if (result instanceof Promise) {
      result.then(done).catch((err) => {
        console.error(`${name} generation failed:`, err)
        process.exitCode = 1
      })
    } else {
      done()
    }
  } catch (error) {
    console.error(`${name} generation failed:`, error)
    process.exitCode = 1
  }
}

/** Uniform "Successfully generated" log with optional parenthesized detail. */
export function logWriteSuccess(targetPath: string, detail?: string): void {
  console.log(`Successfully generated: ${targetPath}${detail ? ` (${detail})` : ''}`)
}

/** Markdown file predicate honoring the -en identity contract per locale. */
export function mdFileFilter(locale: string): (filePath: string) => boolean {
  return (filePath) =>
    locale === 'zh-CN'
      ? /\.md$/i.test(filePath) && !filePath.endsWith('-en.md')
      : /\.md$/i.test(filePath)
}
