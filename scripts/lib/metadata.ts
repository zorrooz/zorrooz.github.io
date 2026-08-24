/**
 * projects/topics 元数据生成骨架：「读 categories.yaml → 归一化 → 写 JSON」。
 * generateProjects/generateTopics 共用，消除生成器互引；配合 lib/yamlEntries.ts。
 */
import fs from 'fs'
import path from 'path'

import { contentDir, localeSuffix, srcDirFor } from '../dataConfig.ts'
import { readYaml, safeArray, writeJsonFile } from './fs.ts'
import { logWriteSuccess } from './cli.ts'
import { normalizeProjectTopicEntry } from './yamlEntries.ts'
import { logError, logWarn } from './log.ts'

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

export type MetadataKind = 'projects' | 'topics'
export type MetadataItem = ProjectItem | TopicItem

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

export function generateMetadataItems(kind: MetadataKind, locale = 'zh-CN') {
  const suffix = localeSuffix(locale)
  /** 源按 locale 分层：zh 读手写源，en 读机器翻译层 */
  const yamlPath = path.join(srcDirFor(locale), `categories${suffix}.yaml`)
  const outputPath = path.join(contentDir, `${kind}${suffix}.json`)

  let items: MetadataItem[] = []
  if (fs.existsSync(yamlPath)) {
    const yamlConfig = readYaml(yamlPath) as Record<string, unknown> | null
    if (yamlConfig === null) {
      logError(`解析 ${kind} 失败: ${yamlPath}:`)
    } else {
      items = safeArray(yamlConfig?.[kind])
        .map((p) => buildMetadataItem(kind, p as Record<string, unknown>))
        .filter((p): p is MetadataItem => p !== null)
    }
  } else {
    logWarn(`categories yaml 不存在: ${yamlPath}`)
  }

  items.sort((a, b) => a.name.localeCompare(b.name))

  try {
    writeJsonFile(outputPath, items)
    logWriteSuccess(outputPath, `${items.length} ${kind}`)
  } catch (e) {
    logError(`写入 ${locale} ${kind}.json 失败:`, e)
  }
}
