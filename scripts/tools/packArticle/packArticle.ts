/**
 * 文章导出核心：把文章 md 及其引用图片导出为「自包含包」。
 * 只读 content-src 源，不修改任何源文件；目标路径每次覆盖重写（快照语义）。
 */
import fs from 'fs'
import path from 'path'
import JSZip from 'jszip'

import { contentSrcDir, dataDir } from '../../dataConfig.ts'

const IMG_RE = /\.(png|jpe?g|gif|svg|webp)([?#].*)?$/i
export const DEFAULT_OUT = path.join(dataDir, 'exports')

export interface ExportResult {
  images: string[]
  outPath: string
}

/**
 * 解析图片引用到 content-src 内的真实文件；外链/非图片返回 null。
 * 支持三种形态：vault 绝对（/assets/x.png）、Obsidian 最短路径（assets/x.png）、
 * 相对路径（../assets/x.png，以文章目录为基准）。
 */
function resolveImageSrc(mdDir: string, href: string): string | null {
  const clean = href.trim()
  if (!IMG_RE.test(clean) || /^(https?:)?\/\//i.test(clean)) return null
  const candidates = clean.startsWith('/')
    ? [path.resolve(contentSrcDir, clean.replace(/^\/+/, ''))]
    : [path.resolve(mdDir, clean), path.resolve(contentSrcDir, clean.replace(/^\.\//, ''))]
  for (const abs of candidates) {
    const rel = path.relative(contentSrcDir, abs)
    if (rel === '..' || rel.startsWith('..')) continue
    if (fs.existsSync(abs) && fs.statSync(abs).isFile()) return abs
  }
  return null
}

/**
 * 导出单篇文章：读源 → 复制引用图 → 改写引用为同目录文件名 → 写导出包（目录或 zip）。
 * 同名图片冲突时第二个起加 slug 前缀；引用中的外链与非图片路径保持原样。
 */
export async function exportArticle(
  rel: string,
  outRoot: string,
  asZip: boolean,
): Promise<ExportResult> {
  const zhPath = path.join(contentSrcDir, `${rel}.md`)
  if (!fs.existsSync(zhPath)) throw new Error(`not found: ${zhPath}`)

  const zhDir = path.dirname(zhPath)
  const slug = path.basename(zhPath, '.md')
  const zhRaw = fs.readFileSync(zhPath, 'utf-8')

  const images = new Map<string, string>()
  const rewritten = zhRaw.replace(
    /!\[([^\]]*)\]\(([^)]+)\)/g,
    (match, alt: string, href: string) => {
      const src = resolveImageSrc(zhDir, href)
      if (!src) return match
      const base = path.basename(src)
      const taken = new Set(images.values())
      const name = taken.has(base) ? `${slug}-${base}` : base
      images.set(src, name)
      return `![${alt}](${name})`
    },
  )

  const relDir = path.dirname(rel)
  if (asZip) {
    const zip = new JSZip()
    zip.file(`${slug}.md`, rewritten)
    for (const [src, name] of images) zip.file(name, fs.readFileSync(src))
    const zipPath = path.join(outRoot, relDir, `${slug}.zip`)
    fs.mkdirSync(path.dirname(zipPath), { recursive: true })
    fs.writeFileSync(zipPath, await zip.generateAsync({ type: 'nodebuffer' }))
    return { images: [...images.keys()], outPath: zipPath }
  }

  const pkgDir = path.join(outRoot, relDir, slug)
  fs.rmSync(pkgDir, { recursive: true, force: true })
  fs.mkdirSync(pkgDir, { recursive: true })
  fs.writeFileSync(path.join(pkgDir, `${slug}.md`), rewritten, 'utf-8')
  for (const [src, name] of images) fs.copyFileSync(src, path.join(pkgDir, name))
  return { images: [...images.keys()], outPath: pkgDir }
}
