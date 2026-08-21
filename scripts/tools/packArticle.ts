// 导出文章：把文章 md 及其引用的图片导出为「自包含包」（默认 ../blog-data/exports/）。
// 与站点部署无关——渲染层原生支持 assets/ 引用，导出是独立的分享/迁移产物。
// 只读 content-src 源，不修改任何源文件；目标路径每次覆盖重写（快照语义）。
// 用法：
//   npm run data:pack -- notes/<cat>/<sub>/<slug> [more...]  # 导出单篇/多篇
//   npm run data:pack                                         # 无参数：导出全部文章
// 选项：
//   --out <dir>   导出根目录（默认 ../blog-data/exports/）
//   --zip         每篇导出为 zip 包（内含 md + 引用图），替代目录包
import fs from 'fs'
import path from 'path'
import { Command } from 'commander'
import JSZip from 'jszip'

import { contentSrcDir, dataDir } from '../dataConfig.ts'
import { mdFileFilter, walk } from '../generators/core/index.ts'

const IMG_RE = /\.(png|jpe?g|gif|svg|webp)([?#].*)?$/i
const DEFAULT_OUT = path.join(dataDir, 'exports')

interface ExportResult {
  images: string[]
  outPath: string
}

/**
 * 解析图片引用到 content-src 内的真实文件；外链/非图片返回 null。
 * 支持三种形态：vault 绝对（/assets/x.png）、Obsidian 最短路径（assets/x.png，以 vault 根为基准）、
 * 相对路径（../assets/x.png 或 ./x.png，以文章目录为基准）。
 */
function resolveImageSrc(mdDir: string, href: string): string | null {
  const clean = href.trim()
  if (!IMG_RE.test(clean) || /^(https?:)?\/\//i.test(clean)) return null
  const candidates = clean.startsWith('/')
    ? [path.resolve(contentSrcDir, clean.replace(/^\/+/, ''))]
    : [
        path.resolve(mdDir, clean),
        path.resolve(contentSrcDir, clean.replace(/^\.\//, '')),
      ]
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
async function exportArticle(rel: string, outRoot: string, asZip: boolean): Promise<ExportResult> {
  const zhPath = path.join(contentSrcDir, `${rel}.md`)
  if (!fs.existsSync(zhPath)) throw new Error(`not found: ${zhPath}`)

  const zhDir = path.dirname(zhPath)
  const slug = path.basename(zhPath, '.md')
  const zhRaw = fs.readFileSync(zhPath, 'utf-8')

  const images = new Map<string, string>() // 源绝对路径 → 导出文件名
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

const program = new Command()
program
  .name('pack')
  .description('导出文章为自包含包（md + 引用图，不改动源文件；默认 exports/，--out 覆盖）')
  .argument(
    '[targets...]',
    '文章相对路径（如 notes/Programming/bash/bash-scripting，可多个）；省略则导出全部文章',
  )
  .option('--out <dir>', '导出根目录（默认 ../blog-data/exports/）')
  .option('--zip', '每篇导出为 zip 包（内含 md + 引用图）')
  .action(async (targets: string[], opts: { out?: string; zip?: boolean }) => {
    let rels: string[]
    if (targets.length > 0) {
      rels = targets.map((t) => t.replace(/\.md$/i, '').replace(/\\/g, '/'))
      for (const rel of rels) {
        if (!rel.startsWith('notes/')) {
          console.error(`[ERROR] 参数必须是 notes/... 相对路径: ${rel}`)
          process.exit(1)
        }
      }
    } else {
      const notesRoot = path.join(contentSrcDir, 'notes')
      rels = walk(notesRoot, mdFileFilter('zh-CN')).map((f) =>
        path.relative(contentSrcDir, f).split(path.sep).join('/').replace(/\.md$/i, ''),
      )
      console.log(`[INFO] 导出全部文章（${rels.length} 篇）`)
    }

    const outRoot = path.resolve(opts.out ?? DEFAULT_OUT)
    const asZip = Boolean(opts.zip)
    let totalImages = 0
    let failed = 0
    for (const rel of rels) {
      try {
        const r = await exportArticle(rel, outRoot, asZip)
        totalImages += r.images.length
        const shown = path.relative(outRoot, r.outPath).split(path.sep).join('/')
        const note = r.images.length > 0 ? `（${r.images.length} 图）` : '（无图）'
        console.log(`[OK] ${rel} → ${shown} ${note}`)
        for (const src of r.images) console.log(`      ${path.relative(contentSrcDir, src)}`)
      } catch (e) {
        failed++
        console.error(`[ERROR] ${rel}: ${(e as Error).message}`)
      }
    }
    console.log(
      `[INFO] 导出完成：${rels.length - failed}/${rels.length} 篇，复制 ${totalImages} 图 → ${outRoot}`,
    )
    if (failed > 0) process.exitCode = 1
  })

program.parse(process.argv)
