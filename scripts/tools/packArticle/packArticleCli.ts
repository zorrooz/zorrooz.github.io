/**
 * data:pack CLI：导出文章为自包含包（md + 引用图，不改动源文件）。
 *
 * 用法：
 *   npm run data:pack -- notes/<cat>/<sub>/<slug> [more...]  # 导出单篇/多篇
 *   npm run data:pack                                        # 无参数：导出全部文章
 * 选项：
 *   --out <dir>   导出根目录（默认 ../blog-data/exports/）
 *   --zip         每篇导出为 zip 包，替代目录包
 */
import path from 'path'
import { Command } from 'commander'

import { contentSrcDir } from '../../dataConfig.ts'
import { mdFileFilter, walk } from '../../lib/fs.ts'
import { DEFAULT_OUT, exportArticle } from './packArticle.ts'

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
