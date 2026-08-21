/**
 * fs/数据/路径域共享实现（自 scripts/generators/core 按职责拆分迁移）。
 * 行为与迁移前逐字节一致；core/index.ts 仅 re-export，消费方 import 路径不变。
 */

import fs from 'fs'
import fsPromises from 'fs/promises'
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

/** 异步递归遍历文件，支持 include/exclude 模式（glob 风格 * 通配，按文件名匹配）；返回绝对路径数组 */
export async function walkAsync(
  dir: string,
  options: { include?: string[]; exclude?: string[]; recursive?: boolean } = {},
): Promise<string[]> {
  // 语义沿袭 translator 原 searchFiles：仅匹配 entry.name（文件名），
  // 目录名同样应用 exclude（如 templates/assets 工具目录不遍历）
  const { include = [], exclude = [], recursive = true } = options
  const files: string[] = []

  async function walkDir(currentDir: string) {
    const entries = await fsPromises.readdir(currentDir, { withFileTypes: true })

    for (const entry of entries) {
      const fullPath = path.join(currentDir, entry.name)

      const isExcluded = exclude.some((pattern) => {
        const regex = new RegExp(pattern.replace('*', '.*').replace('.', '\\.'))
        return regex.test(entry.name)
      })

      if (entry.isDirectory()) {
        if (recursive && !isExcluded) await walkDir(fullPath)
      } else if (!isExcluded) {
        const matchesPattern =
          include.length === 0 ||
          include.some((pattern) => {
            const regex = new RegExp(pattern.replace('*', '.*').replace('.', '\\.'))
            return regex.test(entry.name)
          })
        if (matchesPattern) files.push(fullPath)
      }
    }
  }

  await walkDir(dir)
  return files
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

/** md 文件过滤器：zh 排除 -en 身份后缀，en 全收 */
export function mdFileFilter(locale: string): (filePath: string) => boolean {
  return (filePath) =>
    locale === 'zh-CN'
      ? /\.md$/i.test(filePath) && !filePath.endsWith('-en.md')
      : /\.md$/i.test(filePath)
}

/**
 * 在 categories.json 顶层结构中定位 notes 段。
 * 按结构特征而非标题匹配：标题是翻译产物（zh「笔记」/ en「Notes」），不可作为匹配依据。
 * notes 段特征：items 中至少一项携带非空 categories 数组（子分类）；projects/topics 段为空数组。
 */
export function findNotesSection(categoriesArr: unknown): Record<string, unknown> | null {
  if (!Array.isArray(categoriesArr)) return null
  for (const s of categoriesArr) {
    if (!s || typeof s !== 'object') continue
    const section = s as Record<string, unknown>
    if (!Array.isArray(section.items)) continue
    const hasSubCats = (section.items as unknown[]).some((i) => {
      if (!i || typeof i !== 'object') return false
      const cats = (i as Record<string, unknown>).categories
      return Array.isArray(cats) && cats.length > 0
    })
    if (hasSubCats) return section
  }
  return null
}

/** 遍历 categories.json 的全部文章条目（四层结构 section → item → cat → art） */
export function walkCategoryArticles(
  data: unknown,
  visit: (
    art: Record<string, unknown>,
    ctx: {
      section: Record<string, unknown>
      item: Record<string, unknown>
      cat: Record<string, unknown>
    },
  ) => void,
): void {
  for (const section of Array.isArray(data) ? data : []) {
    if (!section || typeof section !== 'object') continue
    for (const item of safeArray((section as Record<string, unknown>).items)) {
      if (!item || typeof item !== 'object') continue
      for (const cat of safeArray((item as Record<string, unknown>).categories)) {
        if (!cat || typeof cat !== 'object') continue
        for (const art of safeArray((cat as Record<string, unknown>).articles)) {
          if (!art || typeof art !== 'object') continue
          visit(art as Record<string, unknown>, {
            section: section as Record<string, unknown>,
            item: item as Record<string, unknown>,
            cat: cat as Record<string, unknown>,
          })
        }
      }
    }
  }
}

/**
 * 读取并校验 JSON 数组（生成器上游产物前置检查）：
 * 文件缺失 / 解析失败 / 非数组 → 统一报错并置退出码，返回 null。
 */
export function requireJsonArray(
  filePath: string,
  label: string,
  hint: string,
): unknown[] | null {
  if (!fs.existsSync(filePath)) {
    console.error(`${label} not found at ${filePath}. ${hint}`)
    process.exitCode = 1
    return null
  }
  const parsed = readJson(filePath)
  if (parsed === null) {
    console.error(`Failed to read/parse ${label}:`)
    process.exitCode = 1
    return null
  }
  if (!Array.isArray(parsed)) {
    console.error(`${label} is not an array. Abort.`)
    process.exitCode = 1
    return null
  }
  return parsed
}
