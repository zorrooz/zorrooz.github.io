import { TranslationManager, needsTranslation } from './translator.ts'
import { TRANSLATION_CONFIG, TRANSLATION_TARGETS } from './translatorConfig.ts'
import { cacheDir } from '../../dataConfig.ts'
import { Command } from 'commander'
import path from 'path'
import fs from 'fs/promises'
import fsSync from 'fs'

const program = new Command()
const manager = new TranslationManager()

interface TranslationState {
  [sourcePath: string]: string
}

/** 翻译增量状态（键为相对 content-src 的路径，与 translator.ts 写入格式一致） */
function loadState(): TranslationState {
  try {
    return JSON.parse(
      fsSync.readFileSync(path.join(cacheDir, '.translate-state.json'), 'utf-8'),
    ) as TranslationState
  } catch {
    return {}
  }
}

program
  .name('translator')
  .description('模块化文件翻译工具 - 支持MD和YAML文件的批量增量翻译')
  .version('2.0.0')

program
  .command('translate [target]')
  .description('翻译文件或目录（默认：数据分支的 content-src）')
  .option('-f, --force', '重新翻译所有文件（忽略增量检测）')
  .option('-n, --new', '仅翻译新文件（跳过所有已翻译文件）')
  .option('-s, --suffix <suffix>', '翻译文件后缀（默认：-en）')
  .action(async (target = TRANSLATION_TARGETS.CONTENT_SRC, options) => {
    try {
      console.log(`[INFO] 开始翻译: ${target}`)

      const translateOptions = {
        force: options.force || false,
        skipExisting: options.new ? true : !options.force,
        outputSuffix: options.suffix || TRANSLATION_CONFIG.OUTPUT_SUFFIX,
      }

      const results = await manager.translate(target, translateOptions)
      const count = results ? results.length : 0

      if (count === 0) {
        console.log('[INFO] 没有需要翻译的文件（所有文件都已是最新）')
      } else {
        console.log(`[INFO] 翻译完成！处理了 ${count} 个文件`)
      }
    } catch (error) {
      console.error('[ERROR] 翻译失败:', (error as Error).message)
      process.exit(1)
    }
  })

program
  .command('status <target>')
  .description('检查文件翻译状态')
  .action(async (target) => {
    try {
      console.log('[INFO] 检查翻译状态...')

      const stats = await fs.stat(target)
      if (!stats.isDirectory()) {
        console.log('[WARN] 请提供目录路径来检查状态')
        return
      }

      let totalFiles = 0
      let translatedFiles = 0
      let needUpdate = 0

      const state = loadState()

      async function checkDirectory(dir: string) {
        const entries = await fs.readdir(dir, { withFileTypes: true })

        for (const entry of entries) {
          const fullPath = path.join(dir, entry.name)

          if (entry.isDirectory()) {
            await checkDirectory(fullPath)
          } else if (entry.isFile()) {
            const ext = path.extname(entry.name).toLowerCase()
            if (TRANSLATION_CONFIG.SUPPORTED_EXTENSIONS.includes(ext)) {
              totalFiles++

              if (!entry.name.includes(TRANSLATION_CONFIG.OUTPUT_SUFFIX)) {
                const sourcePath = fullPath
                const targetPath = fullPath.replace(
                  ext,
                  `${TRANSLATION_CONFIG.OUTPUT_SUFFIX}${ext}`,
                )

                try {
                  await fs.access(targetPath)
                  translatedFiles++

                  /** 与实际翻译同口径：内容 hash 判定（而非 mtime） */
                  if (await needsTranslation(sourcePath, targetPath, state)) {
                    needUpdate++
                  }
                } catch {
                  /* 目标不存在等访问失败 → 视为需要翻译 */
                }
              }
            }
          }
        }
      }

      await checkDirectory(target)

      console.log('\n[INFO] 翻译状态报告:')
      console.log(`[INFO] 总文件数: ${totalFiles}`)
      console.log(`[INFO] 已翻译文件: ${translatedFiles}`)
      console.log(`[INFO] 需要更新: ${needUpdate}`)
      console.log(`[INFO] 未翻译文件: ${totalFiles - translatedFiles}`)

      if (needUpdate > 0) {
        console.log(`\n[INFO] 建议运行: translator translate "${target}" --force`)
      }
    } catch (error) {
      console.error('[ERROR] 检查状态失败:', (error as Error).message)
      process.exit(1)
    }
  })

program
  .command('default', { isDefault: true })
  .description('默认翻译数据分支的 content-src 目录')
  .action(async () => {
    console.log(`[INFO] 使用默认目录: ${TRANSLATION_TARGETS.CONTENT_SRC}`)
    const manager = new TranslationManager()
    try {
      const results = await manager.translate(TRANSLATION_TARGETS.CONTENT_SRC)
      const count = results ? results.length : 0
      console.log(`[INFO] 翻译完成！成功处理了 ${count} 个文件`)
    } catch (error) {
      console.error('[ERROR] 翻译失败:', (error as Error).message)
      process.exit(1)
    }
  })

program.parse(process.argv)
