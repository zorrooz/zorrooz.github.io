import fs from 'fs'
import path from 'path'
import { contentDir, localeSuffix, srcDirFor } from '../dataConfig.ts'
import { isDirectRun, logWriteSuccess, readYaml, runCliScript, writeJsonFile } from './core/index.ts'

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
              name:
                typeof it?.name === 'string'
                  ? it.name
                  : typeof it?.item === 'string'
                    ? it.item
                    : '',
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
  const suffix = localeSuffix(locale)
  const yamlPath = path.join(srcDirFor(locale), `about${suffix}.yaml`)
  const outputPath = path.join(contentDir, `about${suffix}.json`)

  let raw: unknown = {}
  if (fs.existsSync(yamlPath)) {
    const yamlConfig = readYaml(yamlPath)
    if (yamlConfig !== null && yamlConfig !== undefined) {
      raw = yamlConfig
    } else {
      console.warn(`Warn: failed to parse YAML at ${yamlPath}, using empty object.`)
      raw = {}
    }
  } else {
    console.warn(`Warn: source YAML not found at ${yamlPath}, using empty object.`)
  }

  const result = normalize(raw as Record<string, unknown>)

  try {
    writeJsonFile(outputPath, result)
    logWriteSuccess(outputPath)
  } catch (error) {
    console.error(`Failed to generate ${outputPath}:`, error)
  }
}

if (isDirectRun(import.meta)) {
  runCliScript('about.json', () => {
    generateAboutJson('zh-CN')
    generateAboutJson('en-US')
  })
}

export { generateAboutJson }
