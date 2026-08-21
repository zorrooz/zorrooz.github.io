// 打包文章：把 assets/ 中引用的图片复制到文章 md 同目录，并将引用改写为同目录相对路径。
// 用法：
//   npm run data:pack -- notes/<cat>/<sub>/<slug>   # 打包单篇（可传多个）
//   npm run data:pack                                # 无参数：递归打包全部文章（跳过 -en 遗留）
// 幂等：已打包（同目录引用）的文章再次运行不会重复复制。
import fs from 'fs'
import path from 'path'
import { Command } from 'commander'

import { contentSrcDir, enSrcDir } from './dataConfig.ts'
import { mdFileFilter, walk } from './generators/core/index.ts'

const IMG_RE = /\.(png|jpe?g|gif|svg|webp)([?#].*)?$/i

interface PackResult {
  copied: string[]
  rewrote: number
}

/** 打包单篇文章：assets 引用图 → 文章同目录 + 改写引用（同步改写 cache/en 镜像） */
function packArticle(rel: string): PackResult {
  const zhPath = path.join(contentSrcDir, `${rel}.md`)
  if (!fs.existsSync(zhPath)) throw new Error(`not found: ${zhPath}`)

  const zhRaw = fs.readFileSync(zhPath, 'utf-8')
  const zhDir = path.dirname(zhPath)
  const slug = path.basename(zhPath, '.md')
  const mapping = new Map<string, string>() // 原始 href → 同目录新文件名
  const copied: string[] = []

  const rewritten = zhRaw.replace(
    /!\[([^\]]*)\]\(([^)]+)\)/g,
    (match, alt: string, href: string) => {
      const src = resolveAssetSrc(zhDir, href.trim())
      if (!src) return match
      const base = path.basename(src)
      let target = path.join(zhDir, base)
      const sameFile = fs.existsSync(target) && fs.realpathSync(target) === fs.realpathSync(src)
      if (!sameFile && fs.existsSync(target)) {
        // 同目录已存在不同文件：加 slug 前缀避免跨文章冲突
        target = path.join(zhDir, `${slug}-${base}`)
      }
      if (!sameFile) {
        fs.copyFileSync(src, target)
        copied.push(target)
      }
      const newName = path.basename(target)
      mapping.set(href.trim(), newName)
      return `![${alt}](${newName})`
    },
  )

  if (rewritten !== zhRaw) fs.writeFileSync(zhPath, rewritten, 'utf-8')

  // en 镜像（cache/en）：仅改写引用为同目录路径，不复制文件
  //（渲染时 rewriteImageLinks 会解析到 content-src 同目录的图）
  const enPath = path.join(enSrcDir, `${rel}-en.md`)
  if (fs.existsSync(enPath)) {
    const enRaw = fs.readFileSync(enPath, 'utf-8')
    const enRewritten = enRaw.replace(
      /!\[([^\]]*)\]\(([^)]+)\)/g,
      (match, alt: string, href: string) => {
        const newName = mapping.get(href.trim())
        return newName ? `![${alt}](${newName})` : match
      },
    )
    if (enRewritten !== enRaw) fs.writeFileSync(enPath, enRewritten, 'utf-8')
  }

  return { copied, rewrote: mapping.size }
}

/** 解析图片引用到 content-src/assets 下的真实文件；非 assets 引用返回 null */
function resolveAssetSrc(mdDir: string, href: string): string | null {
  if (!IMG_RE.test(href) || /^(https?:)?\/\//i.test(href)) return null
  const abs = href.startsWith('/')
    ? path.resolve(contentSrcDir, href.replace(/^\/+/, ''))
    : path.resolve(mdDir, href)
  const rel = path.relative(contentSrcDir, abs)
  if (rel === '..' || rel.startsWith('..') || !rel.startsWith('assets/')) return null
  return fs.existsSync(abs) ? abs : null
}

const program = new Command()
program
  .name('pack')
  .description('打包文章：把 assets/ 引用图复制到文章同目录并改写引用（幂等）')
  .argument(
    '[targets...]',
    '文章相对路径（如 notes/Programming/bash/bash-scripting，可多个）；省略则打包全部文章',
  )
  .action((targets: string[]) => {
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
      console.log(`[INFO] 递归打包全部文章（${rels.length} 篇）`)
    }

    let totalCopied = 0
    let totalRewrote = 0
    let failed = 0
    for (const rel of rels) {
      try {
        const r = packArticle(rel)
        totalCopied += r.copied.length
        totalRewrote += r.rewrote
        if (r.copied.length === 0 && r.rewrote === 0) {
          console.log(`[INFO] ${rel}（无 assets 图片引用）`)
        } else {
          console.log(`[OK] ${rel}: 复制 ${r.copied.length} 图 / 改写 ${r.rewrote} 引用`)
          for (const f of r.copied) console.log(`      ${path.relative(contentSrcDir, f)}`)
        }
      } catch (e) {
        failed++
        console.error(`[ERROR] ${rel}: ${(e as Error).message}`)
      }
    }
    console.log(
      `[INFO] 打包完成：${rels.length - failed}/${rels.length} 篇，复制 ${totalCopied} 图，改写 ${totalRewrote} 引用`,
    )
    if (failed > 0) process.exitCode = 1
  })

program.parse(process.argv)
