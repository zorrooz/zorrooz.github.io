import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'
import yaml from 'js-yaml'

import { contentSrcDir, contentDir } from '../dataConfig.ts'
import { ensureDirectoryExistence } from './core/index.ts'

/**
 * 项目已改为纯元信息模式：不再从 md 扫描文章，
 * 仅从 categories.yaml 读取项目元信息（md 源文件保留但不生成文章）。
 */
export interface ProjectItem {
  name: string
  title: string
  desc: string
  date: string
  tags: string[]
  github: string
  url: string
  status: string
  language: string
  stars: number
  license: string
  version: string
}

function getFilePaths(locale = 'zh-CN') {
  const suffix = locale === 'zh-CN' ? '' : '-en'
  return {
    yamlPath: path.join(contentSrcDir, `categories${suffix}.yaml`),
    outputPath: path.join(contentDir, `projects${suffix}.json`),
  }
}

function safeArray(x: unknown): unknown[] {
  return Array.isArray(x) ? x : []
}

function buildProjectItem(rawDef: Record<string, unknown>): ProjectItem | null {
  if (typeof rawDef?.name !== 'string' || !rawDef.name) return null
  return {
    name: rawDef.name,
    title: typeof rawDef.title === 'string' ? rawDef.title : rawDef.name,
    desc: typeof rawDef.desc === 'string' ? rawDef.desc : '',
    date: typeof rawDef.date === 'string' ? rawDef.date : '',
    tags: safeArray(rawDef.tags).filter((t) => typeof t === 'string'),
    github: typeof rawDef.github === 'string' ? rawDef.github : '',
    url: typeof rawDef.url === 'string' ? rawDef.url : '',
    status: typeof rawDef.status === 'string' ? rawDef.status : '',
    language: typeof rawDef.language === 'string' ? rawDef.language : '',
    stars: Number(rawDef.stars) || 0,
    license: typeof rawDef.license === 'string' ? rawDef.license : '',
    version: typeof rawDef.version === 'string' ? rawDef.version : '',
  }
}

function generateProjectsJson(locale = 'zh-CN') {
  const paths = getFilePaths(locale)

  let items: ProjectItem[] = []
  try {
    if (fs.existsSync(paths.yamlPath)) {
      const raw = fs.readFileSync(paths.yamlPath, 'utf-8')
      const obj = yaml.load(raw) as Record<string, unknown> | undefined
      items = safeArray(obj?.projects)
        .map((p) => buildProjectItem(p as Record<string, unknown>))
        .filter((p): p is ProjectItem => p !== null)
    } else {
      console.warn(`Warn: categories yaml not found at ${paths.yamlPath}`)
    }
  } catch (e) {
    console.error(`Failed to parse projects from ${paths.yamlPath}:`, e)
  }

  items.sort((a, b) => a.name.localeCompare(b.name))

  try {
    const targetPath = paths.outputPath
    ensureDirectoryExistence(targetPath)
    fs.writeFileSync(targetPath, JSON.stringify(items, null, 2), 'utf-8')
    console.log(`Successfully generated: ${targetPath} (${items.length} projects)`)
  } catch (e) {
    console.error(`Failed to write ${locale} projects.json:`, e)
  }
}

function main() {
  console.log('Starting projects.json generation script...')
  generateProjectsJson('zh-CN')
  generateProjectsJson('en-US')
  console.log('projects.json generation complete.')
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main()
}

export { generateProjectsJson }
