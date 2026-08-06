import fs from 'fs'
import path from 'path'

import { contentDir, localeSuffix, srcDirFor } from '../dataConfig.ts'
import {
  readYaml,
  safeArray,
  writeJsonFile,
  isDirectRun,
  runCliScript,
  logWriteSuccess,
} from './core/index.ts'

/**
 * 项目为纯元信息模式：不从 md 扫描文章，仅从 categories.yaml 读取项目元数据。
 */
interface ProjectItem {
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
  const suffix = localeSuffix(locale)
  return {
    // 源按 locale 分层：zh 读手写源，en 读机器翻译层
    yamlPath: path.join(srcDirFor(locale), `categories${suffix}.yaml`),
    outputPath: path.join(contentDir, `projects${suffix}.json`),
  }
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
  if (fs.existsSync(paths.yamlPath)) {
    const yamlConfig = readYaml(paths.yamlPath) as Record<string, unknown> | null
    if (yamlConfig === null) {
      console.error(`Failed to parse projects from ${paths.yamlPath}:`)
    } else {
      items = safeArray(yamlConfig?.projects)
        .map((p) => buildProjectItem(p as Record<string, unknown>))
        .filter((p): p is ProjectItem => p !== null)
    }
  } else {
    console.warn(`Warn: categories yaml not found at ${paths.yamlPath}`)
  }

  items.sort((a, b) => a.name.localeCompare(b.name))

  try {
    const targetPath = paths.outputPath
    writeJsonFile(targetPath, items)
    logWriteSuccess(targetPath, `${items.length} projects`)
  } catch (e) {
    console.error(`Failed to write ${locale} projects.json:`, e)
  }
}

if (isDirectRun(import.meta)) {
  runCliScript('projects.json', () => {
    generateProjectsJson('zh-CN')
    generateProjectsJson('en-US')
  })
}

export { generateProjectsJson }
