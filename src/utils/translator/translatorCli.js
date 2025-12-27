import { TranslationManager } from './translator.js'
import { TRANSLATION_CONFIG } from './translatorConfig.js'
import { Command } from 'commander'
import path from 'path'
import fs from 'fs/promises'

const program = new Command()
const manager = new TranslationManager()

// 设置命令行程序
program
  .name('translator')
  .description('模块化文件翻译工具 - 支持MD和YAML文件的批量增量翻译')
  .version('2.0.0')

// 翻译命令
program
  .command('translate [target]')
  .description('翻译文件或目录（默认：src/content-src）')
  .option('-f, --force', '重新翻译所有文件（忽略增量检测）')
  .option('-n, --new', '仅翻译新文件（跳过所有已翻译文件）')
  .option('-s, --suffix <suffix>', '翻译文件后缀（默认：-en）')
  .action(async (target = 'src/content-src', options) => {
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
      console.error('[ERROR] 翻译失败:', error.message)
      process.exit(1)
    }
  })

// 检查翻译状态
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

      async function checkDirectory(dir) {
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

                  const sourceStats = await fs.stat(sourcePath)
                  const targetStats = await fs.stat(targetPath)
                  if (sourceStats.mtime > targetStats.mtime) {
                    needUpdate++
                  }
                } catch {
                  // 翻译文件不存在
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
      console.error('[ERROR] 检查状态失败:', error.message)
      process.exit(1)
    }
  })

// 显示帮助信息
program
  .command('help')
  .description('显示帮助信息')
  .action(() => {
    console.log(`
命令:
  translate [目录路径]    翻译文件或目录
  status <目录路径>       检查翻译状态
  help                    显示此帮助信息

选项:
  -f, --force            强制重新翻译所有文件
  -n, --new              仅翻译新文件
  -s, --suffix <后缀>    翻译文件后缀（默认：-en）

示例:
  npm run translate                                     # 增量翻译默认目录
  npm run translate -- translate                        # 同上
  npm run translate -- translate --force                # 强制重新翻译
  npm run translate -- translate src/content-src --new  # 仅翻译新文件
  npm run translate -- status src/content-src           # 检查状态
    `)
  })

// 处理未知命令
program.on('command:*', () => {
  console.error('[ERROR] 未知命令: %s', program.args.join(' '))
  console.log('请使用 "npm run translate -- help" 查看可用命令')
  process.exit(1)
})

// 默认命令 - 当没有提供子命令时执行
program
  .command('default', { isDefault: true })
  .description('默认翻译 src/content-src 目录')
  .action(async () => {
    console.log('[INFO] 使用默认目录: src/content-src')
    const manager = new TranslationManager()
    try {
      const results = await manager.translate('src/content-src')
      const count = results ? results.length : 0
      console.log(`[INFO] 翻译完成！成功处理了 ${count} 个文件`)
    } catch (error) {
      console.error('[ERROR] 翻译失败:', error.message)
      process.exit(1)
    }
  })

// 解析命令行参数
program.parse(process.argv)
