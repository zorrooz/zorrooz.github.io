import fs from 'fs'
import path from 'path'

import { contentDir, localeSuffix, srcDirFor } from '../dataConfig.ts'
import { normalizeProjectTopicEntry } from '../lib/yamlEntries.ts'
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
 * projects/topics 共用同一套「读 categories.yaml → 归一化 → 写 JSON」骨架，
 * 仅条目字段集不同，按 kind 分发。
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

interface TopicItem {
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

type MetadataKind = 'projects' | 'topics'
type MetadataItem = ProjectItem | TopicItem

function buildMetadataItem(
  kind: MetadataKind,
  rawDef: Record<string, unknown>,
): MetadataItem | null {
  const item = normalizeProjectTopicEntry(rawDef)
  if (!item.name) return null
  if (kind === 'projects') {
    return {
      name: item.name,
      title: item.title,
      desc: item.desc,
      date: item.date,
      tags: item.tags,
      github: item.github,
      url: item.url,
      status: item.status,
      language: item.language,
      stars: item.stars,
      license: item.license,
      version: item.version,
    }
  }
  return {
    name: item.name,
    title: item.title,
    desc: item.desc,
    date: item.date,
    tags: item.tags,
    doi: item.doi,
    url: item.url,
    status: item.status,
    journal: item.journal,
    year: item.year,
    authors: item.authors,
  }
}

function generateMetadataItems(kind: MetadataKind, locale = 'zh-CN') {
  const suffix = localeSuffix(locale)
  // 源按 locale 分层：zh 读手写源，en 读机器翻译层
  const yamlPath = path.join(srcDirFor(locale), `categories${suffix}.yaml`)
  const outputPath = path.join(contentDir, `${kind}${suffix}.json`)

  let items: MetadataItem[] = []
  if (fs.existsSync(yamlPath)) {
    const yamlConfig = readYaml(yamlPath) as Record<string, unknown> | null
    if (yamlConfig === null) {
      console.error(`Failed to parse ${kind} from ${yamlPath}:`)
    } else {
      items = safeArray(yamlConfig?.[kind])
        .map((p) => buildMetadataItem(kind, p as Record<string, unknown>))
        .filter((p): p is MetadataItem => p !== null)
    }
  } else {
    console.warn(`Warn: categories yaml not found at ${yamlPath}`)
  }

  items.sort((a, b) => a.name.localeCompare(b.name))

  try {
    writeJsonFile(outputPath, items)
    logWriteSuccess(outputPath, `${items.length} ${kind}`)
  } catch (e) {
    console.error(`Failed to write ${locale} ${kind}.json:`, e)
  }
}

function generateProjectsJson(locale = 'zh-CN') {
  generateMetadataItems('projects', locale)
}

if (isDirectRun(import.meta)) {
  runCliScript('projects.json', () => {
    generateProjectsJson('zh-CN')
    generateProjectsJson('en-US')
  })
}

export { generateMetadataItems, generateProjectsJson }
