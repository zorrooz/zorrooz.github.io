/**
 * TagMerger CLI — AI 辅助标签合并工具
 *
 * 用法：
 *   node ./src/utils/tagMerger/tagMergerCli.ts                # dry-run：生成映射建议
 *   node ./src/utils/tagMerger/tagMergerCli.ts --apply        # 将映射写回源文件
 *   node ./src/utils/tagMerger/tagMergerCli.ts --locale en    # 处理英文标签
 *   node ./src/utils/tagMerger/tagMergerCli.ts --force        # 忽略已有映射文件重新生成
 */

import fs from 'fs'
import path from 'path'
import { program } from 'commander'

import {
  collectTags,
  formatTagStats,
  generateMappingWithLlm,
  loadMapping,
  applyMappingToMarkdown,
  applyMappingToCategoriesYaml,
  defaultMappingPath,
} from './tagMerger.ts'
import { loadLlmConfig } from './tagMerger.ts'

program
  .name('tagmerge')
  .description('AI 辅助标签合并工具：收集全部标签，LLM 生成合并映射，可写回源文件')
  .option('--dry-run', '仅生成映射建议文件，不写回源文件（默认）')
  .option('--apply', '生成映射并写回 md frontmatter 与 categories.yaml')
  .option('--locale <lang>', '处理语言：zh（默认）或 en', 'zh')
  .option('--output <path>', `映射文件输出路径（默认 ${defaultMappingPath}）`, defaultMappingPath)
  .option('--force', '忽略已存在的映射文件，强制重新生成')
  .parse(process.argv)

const opts = program.opts()

function resolveLocale(raw: string): 'zh-CN' | 'en-US' {
  const r = raw.toLowerCase()
  if (r === 'zh' || r === 'zh-cn') return 'zh-CN'
  if (r === 'en' || r === 'en-us') return 'en-US'
  console.error('[ERROR] 无效的 --locale 值，仅支持 zh / en（或 zh-CN / en-US）')
  process.exit(1)
}

const locale = resolveLocale(String(opts.locale || 'zh'))
const outputPath = path.resolve(opts.output)
const apply = Boolean(opts.apply)

async function main() {
  console.log(`== TagMerger (${locale}) ==`)

  const tags = collectTags(locale)
  if (tags.size === 0) {
    console.error('[ERROR] 未收集到任何标签')
    process.exitCode = 1
    return
  }
  console.log(`收集到 ${tags.size} 个标签：`)
  console.log(formatTagStats(tags))

  let mapping: Record<string, string>
  const existing = loadMapping(outputPath)
  if (existing && existing.locale === locale && !opts.force) {
    console.log(`[INFO] 使用已有映射文件：${outputPath}（--force 可重新生成）`)
    mapping = existing.mapping
  } else {
    const llm = await loadLlmConfig()
    if (!llm) {
      console.error(
        '[ERROR] 未找到 LLM 配置：请在 src/config/llmConfig.ts 配置 { url, apikey, model }，该文件被 gitignore。',
      )
      process.exitCode = 1
      return
    }
    console.log('[INFO] 正在调用 LLM 生成合并映射...')
    mapping = await generateMappingWithLlm(tags, locale, llm)
    const out: { locale: string; generatedAt: string; mapping: Record<string, string> } = {
      locale,
      generatedAt: new Date().toISOString(),
      mapping,
    }
    fs.writeFileSync(outputPath, JSON.stringify(out, null, 2), 'utf-8')
    console.log(`[INFO] 映射文件已写入：${outputPath}`)
  }

  const changes = Object.entries(mapping).filter(([oldTag, newTag]) => oldTag !== newTag)
  console.log(`\n合并建议（${changes.length} 项）：`)
  for (const [oldTag, newTag] of changes) {
    console.log(`  ${oldTag} -> ${newTag}`)
  }
  if (changes.length === 0) {
    console.log('（无合并项，所有标签已是规范形式）')
  }

  if (apply) {
    const mdChanged = applyMappingToMarkdown(mapping, locale)
    const yamlChanged = applyMappingToCategoriesYaml(mapping, locale)
    console.log(`[INFO] 已写回：${mdChanged} 个 md 文件、${yamlChanged} 个 categories.yaml`)
  } else {
    console.log(`[INFO] dry-run 模式，未写回源文件。确认映射后使用 --apply 应用。`)
  }
}

main().catch((err) => {
  console.error('[ERROR]', err instanceof Error ? err.message : err)
  process.exitCode = 1
})
