import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'
import yaml from 'js-yaml'
import { contentSrcDir, contentDir } from '../dataConfig.ts'
import { ensureDirectoryExistence } from './core/index.ts'

const isDirectRun =
  process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)

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
    ? raw.items.map((it: Record<string, unknown>): ResourceItem => ({
        name: typeof it?.name === 'string' ? it.name : '',
        url: typeof it?.url === 'string' ? it.url : '',
        desc: typeof it?.desc === 'string' ? it.desc : '',
      }))
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
  const yamlPath = path.join(
    contentSrcDir,
    locale === 'zh-CN' ? 'resources.yaml' : 'resources-en.yaml',
  )
  const outputPath = path.join(
    contentDir,
    locale === 'zh-CN' ? 'resources.json' : 'resources-en.json',
  )

  let raw: unknown = []
  if (fs.existsSync(yamlPath)) {
    try {
      const yamlContent = fs.readFileSync(yamlPath, 'utf-8')
      const parsed = yaml.load(yamlContent)
      raw = parsed || []
    } catch (e) {
      console.warn(
        `Warn: failed to parse YAML at ${yamlPath}, using empty array.`,
        e instanceof Error ? e.message : e,
      )
      raw = []
    }
  } else {
    console.warn(`Warn: source YAML not found at ${yamlPath}, using empty array.`)
  }

  const result = normalize(raw)

  try {
    ensureDirectoryExistence(outputPath)
    fs.writeFileSync(outputPath, JSON.stringify(result, null, 2), 'utf-8')
    console.log(`Successfully generated: ${outputPath}`)
  } catch (error) {
    console.error(`Failed to generate ${outputPath}:`, error)
  }
}

function main() {
  console.log('Starting resources.json generation script...')
  generateResourcesJson('zh-CN')
  generateResourcesJson('en-US')
  console.log('resources.json generation complete.')
}

if (isDirectRun) main()

export { generateResourcesJson }
