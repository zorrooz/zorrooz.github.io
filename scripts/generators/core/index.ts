/**
 * 内容生成器共享工具。行为须与重构前各生成器的原实现逐字节一致，
 * 重构后须对比前后生成的 JSON 验证。
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

/** 粗略去除 markdown 标记 → 纯文本，用于字数统计 */
export function markdownToPlain(text: string): string {
  // 先剥代码块与行内代码
  let t = text.replace(/```[\s\S]*?```/g, ' ')
  t = t.replace(/`[^`]*`/g, ' ')
  // 图像/链接/HTML 标签
  t = t.replace(/!\[[^\]]*]\([^)]+\)/g, ' ')
  t = t.replace(/\[[^\]]*]\([^)]+\)/g, ' ')
  t = t.replace(/<[^>]+>/g, ' ')
  // 标题、加粗/斜体/引用/列表/表格标记与参考链接
  t = t.replace(/^#{1,6}\s+/gm, ' ')
  t = t.replace(/[*_~`>#|-]{1,}/g, ' ')
  t = t.replace(/^\s*\[[^\]]+]:\s+\S+.*$/gm, ' ')
  // 折叠空白
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
 * 读取 JSON 文件：文件缺失或非法返回 null（文件内容为字面 `null` 时返回 null 本身）。
 * 调用方可用 fs.existsSync 区分「缺失」。
 */
export function readJson(filePath: string): unknown {
  try {
    if (!fs.existsSync(filePath)) return null
    return JSON.parse(fs.readFileSync(filePath, 'utf-8').trim())
  } catch {
    return null
  }
}

/** 读取 YAML 文件；缺失返回 null，空文件返回 undefined */
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
 * 判断模块是否被直接执行（node foo.ts），容忍相对/绝对 argv 形式。
 * 用于门控 CLI 入口，避免 import 生成器模块时误触发其 main()。
 */
export function isDirectRun(importMeta: ImportMeta): boolean {
  return Boolean(process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(importMeta.url))
}

/** 单个生成器 CLI 的通用引导：banner + 执行 + 失败退出码 */
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

/** 统一的「生成成功」日志，可附括号明细 */
export function logWriteSuccess(targetPath: string, detail?: string): void {
  console.log(`Successfully generated: ${targetPath}${detail ? ` (${detail})` : ''}`)
}

/** md 文件过滤器：zh 排除 -en 身份后缀，en 全收 */
export function mdFileFilter(locale: string): (filePath: string) => boolean {
  return (filePath) =>
    locale === 'zh-CN'
      ? /\.md$/i.test(filePath) && !filePath.endsWith('-en.md')
      : /\.md$/i.test(filePath)
}
