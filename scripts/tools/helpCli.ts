/**
 * 统一工具帮助入口：npm run help [<name>]
 * 无参列出全部数据工具；带名字则转发到对应 CLI 的 --help。
 */
import { spawnSync } from 'child_process'
import path from 'path'
import { fileURLToPath } from 'url'

const here = path.dirname(fileURLToPath(import.meta.url))

interface ToolEntry {
  desc: string
  /** 相对于本目录的 CLI 路径；缺省 = 非 CLI（git 命令串），仅展示说明 */
  cli?: string
}

const TOOLS: Record<string, ToolEntry> = {
  translate: { desc: 'AI 增量翻译（content-src → cache/en）', cli: 'translator/translatorCli.ts' },
  tagmerge: { desc: 'zh→en 标签映射增量补齐', cli: 'tagMerger/tagMergerCli.ts' },
  pack: { desc: '导出文章自包含包（md + 引用图）', cli: 'packArticle/packArticleCli.ts' },
  deploy: { desc: 'git add+commit+push 数据分支（触发 CI 部署）' },
  publish: { desc: '仅 git push 数据分支（触发 CI 部署）' },
}

function listAll() {
  console.log('数据工具一览：')
  for (const [name, { desc }] of Object.entries(TOOLS)) {
    console.log(`  ${name.padEnd(10)} ${desc}`)
  }
  console.log('\n用法：npm run help <name>   例如 npm run help translate')
}

const name = process.argv[2]
if (!name) {
  listAll()
} else {
  const entry = TOOLS[name]
  if (!entry) {
    console.error(`未知工具：${name}`)
    listAll()
    process.exitCode = 1
  } else if (entry.cli) {
    const result = spawnSync(process.execPath, [path.resolve(here, entry.cli), '--help'], {
      stdio: 'inherit',
    })
    if (result.error || (result.status ?? 1) !== 0) process.exitCode = result.status ?? 1
  } else {
    console.log(`${name} — ${entry.desc}\n（git 命令串，无子命令帮助；见 package.json scripts 与 AGENTS.md）`)
  }
}
