import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'
import yaml from 'js-yaml'

import { contentDir, localeSuffix, srcDirFor } from '../dataConfig.ts'
import { ensureDirectoryExistence } from './core/index.ts'

/**
 * 课题已改为纯元信息模式：不再从 md 扫描文章，
 * 仅从 categories.yaml 读取课题元信息（md 源文件保留但不生成文章）。
 */
export interface TopicItem {
  name: string
  title: string
  desc: string
  date: string
  tags: string[]
  doi: string
  url: string
  status: string
  journal: string
  year: number
  authors: string[]
}

function getFilePaths(locale = 'zh-CN') {
  const suffix = localeSuffix(locale)
  return {
    // 源按 locale 分层：zh 读手写源，en 读机器翻译层
    yamlPath: path.join(srcDirFor(locale), `categories${suffix}.yaml`),
    outputPath: path.join(contentDir, `topics${suffix}.json`),
  }
}

function safeArray(x: unknown): unknown[] {
  return Array.isArray(x) ? x : []
}

function buildTopicItem(rawDef: Record<string, unknown>): TopicItem | null {
  if (typeof rawDef?.name !== 'string' || !rawDef.name) return null
  return {
    name: rawDef.name,
    title: typeof rawDef.title === 'string' ? rawDef.title : rawDef.name,
    desc: typeof rawDef.desc === 'string' ? rawDef.desc : '',
    date: typeof rawDef.date === 'string' ? rawDef.date : '',
    tags: safeArray(rawDef.tags).filter((t) => typeof t === 'string'),
    doi: typeof rawDef.doi === 'string' ? rawDef.doi : '',
    url: typeof rawDef.url === 'string' ? rawDef.url : '',
    status: typeof rawDef.status === 'string' ? rawDef.status : '',
    journal: typeof rawDef.journal === 'string' ? rawDef.journal : '',
    year: Number(rawDef.year) || 0,
    authors: safeArray(rawDef.authors).filter((t) => typeof t === 'string'),
  }
}

function generateTopicsJson(locale = 'zh-CN') {
  const paths = getFilePaths(locale)

  let items: TopicItem[] = []
  try {
    if (fs.existsSync(paths.yamlPath)) {
      const raw = fs.readFileSync(paths.yamlPath, 'utf-8')
      const obj = yaml.load(raw) as Record<string, unknown> | undefined
      items = safeArray(obj?.topics)
        .map((t) => buildTopicItem(t as Record<string, unknown>))
        .filter((t): t is TopicItem => t !== null)
    } else {
      console.warn(`Warn: categories yaml not found at ${paths.yamlPath}`)
    }
  } catch (e) {
    console.error(`Failed to parse topics from ${paths.yamlPath}:`, e)
  }

  items.sort((a, b) => a.name.localeCompare(b.name))

  try {
    const targetPath = paths.outputPath
    ensureDirectoryExistence(targetPath)
    fs.writeFileSync(targetPath, JSON.stringify(items, null, 2), 'utf-8')
    console.log(`Successfully generated: ${targetPath} (${items.length} topics)`)
  } catch (e) {
    console.error(`Failed to write ${locale} topics.json:`, e)
  }
}

function main() {
  console.log('Starting topics.json generation script...')
  generateTopicsJson('zh-CN')
  generateTopicsJson('en-US')
  console.log('topics.json generation complete.')
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main()
}

export { generateTopicsJson }
