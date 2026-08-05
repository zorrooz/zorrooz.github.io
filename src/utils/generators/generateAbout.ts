import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'
import yaml from 'js-yaml'
import { contentSrcDir, contentDir } from '../dataConfig.ts'
import { ensureDirectoryExistence } from './core/index.ts'

const isDirectRun =
  process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)

interface ExperienceItem {
  year: string
  title: string
  desc: string
}

interface AboutSectionItem {
  name: string
  desc: string
}

interface AboutSection {
  title: string
  items: AboutSectionItem[]
}

interface AboutData {
  introduction: string
  experience: ExperienceItem[]
  section: AboutSection[]
  contacts: unknown[]
}

function normalize(raw: Record<string, unknown> = {}): AboutData {
  const intro = typeof raw.introduction === 'string' ? raw.introduction : ''

  const experience: ExperienceItem[] = Array.isArray(raw.experience)
    ? raw.experience.map((it: Record<string, unknown>) => ({
        year: typeof it?.year === 'string' ? it.year : '',
        title: typeof it?.title === 'string' ? it.title : '',
        desc: typeof it?.desc === 'string' ? it.desc : '',
      }))
    : []

  const section: AboutSection[] = Array.isArray(raw.section)
    ? raw.section.map((s: Record<string, unknown>) => {
        const title = typeof s?.title === 'string' ? s.title : ''
        const items = Array.isArray(s?.items)
          ? s.items.map((it: Record<string, unknown>) => ({
              name: typeof it?.name === 'string' ? it.name : typeof it?.item === 'string' ? it.item : '',
              desc: typeof it?.desc === 'string' ? it.desc : '',
            }))
          : []
        return { title, items }
      })
    : []

  const contacts = Array.isArray(raw.contacts) ? raw.contacts : []

  return {
    introduction: intro,
    experience,
    section,
    contacts,
  }
}

function generateAboutJson(locale = 'zh-CN') {
  const yamlPath = path.join(contentSrcDir, locale === 'zh-CN' ? 'about.yaml' : 'about-en.yaml')
  const outputPath = path.join(
    contentDir,
    locale === 'zh-CN' ? 'about.json' : 'about-en.json',
  )

  let raw: Record<string, unknown> = {}
  if (fs.existsSync(yamlPath)) {
    try {
      const yamlContent = fs.readFileSync(yamlPath, 'utf-8')
      raw = (yaml.load(yamlContent) as Record<string, unknown>) || {}
    } catch (e) {
      console.warn(
        `Warn: failed to parse YAML at ${yamlPath}, using empty object.`,
        e instanceof Error ? e.message : e,
      )
      raw = {}
    }
  } else {
    console.warn(`Warn: source YAML not found at ${yamlPath}, using empty object.`)
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
  console.log('Starting about.json generation script...')
  generateAboutJson('zh-CN')
  generateAboutJson('en-US')
  console.log('about.json generation complete.')
}

if (isDirectRun) main()

export { generateAboutJson }
