import fs from 'fs'
import path from 'path'
import { contentDir, localeSuffix, srcDirFor } from '../dataConfig.ts'
import {
  isDirectRun,
  logWriteSuccess,
  readYaml,
  runCliScript,
  writeJsonFile,
} from './core/index.ts'

interface ResourceItem {
  name: string
  url: string
  desc: string
}

interface ResourceNode {
  title: string
  items?: ResourceItem[]
  children?: ResourceNode[]
}

interface ResourceCategory {
  title: string
  children: ResourceNode[]
}

function normalizeNode(raw: Record<string, unknown>): ResourceNode {
  const title = typeof raw?.title === 'string' ? raw.title : ''
  const items = Array.isArray(raw?.items)
    ? raw.items.map(
        (it: Record<string, unknown>): ResourceItem => ({
          name: typeof it?.name === 'string' ? it.name : '',
          url: typeof it?.url === 'string' ? it.url : '',
          desc: typeof it?.desc === 'string' ? it.desc : '',
        }),
      )
    : []
  const children = Array.isArray(raw?.children) ? raw.children.map(normalizeNode) : []
  const node: ResourceNode = { title }
  if (items.length) node.items = items
  if (children.length) node.children = children
  return node
}

function normalize(raw: unknown): ResourceCategory[] {
  const list = Array.isArray(raw) ? raw : []
  return list.map((cat: Record<string, unknown>) => {
    const title = typeof cat?.title === 'string' ? cat.title : ''
    const children = Array.isArray(cat?.children) ? cat.children : []
    return { title, children: children.map(normalizeNode) }
  })
}

function generateResourcesJson(locale = 'zh-CN') {
  const suffix = localeSuffix(locale)
  const yamlPath = path.join(srcDirFor(locale), `resources${suffix}.yaml`)
  const outputPath = path.join(contentDir, `resources${suffix}.json`)

  let raw: unknown = []
  if (fs.existsSync(yamlPath)) {
    const yamlConfig = readYaml(yamlPath)
    if (yamlConfig !== null && yamlConfig !== undefined) {
      raw = yamlConfig
    } else {
      console.warn(`Warn: failed to parse YAML at ${yamlPath}, using empty array.`)
      raw = []
    }
  } else {
    console.warn(`Warn: source YAML not found at ${yamlPath}, using empty array.`)
  }

  const result = normalize(raw)

  try {
    writeJsonFile(outputPath, result)
    logWriteSuccess(outputPath)
  } catch (error) {
    console.error(`Failed to generate ${outputPath}:`, error)
  }
}

if (isDirectRun(import.meta)) {
  runCliScript('resources.json', () => {
    generateResourcesJson('zh-CN')
    generateResourcesJson('en-US')
  })
}

export { generateResourcesJson }
