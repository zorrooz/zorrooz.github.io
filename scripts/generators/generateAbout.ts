import { isDirectRun, runCliScript } from './core/index.ts'
import { generateJsonFromYaml } from '../lib/yamlJson.ts'
import type {
  AboutContact,
  AboutData,
  AboutExperience,
  AboutSection,
  AboutSectionItem,
} from '../../src/types.ts'

function normalize(raw: Record<string, unknown> = {}): AboutData {
  const intro = typeof raw.introduction === 'string' ? raw.introduction : ''

  const experience: AboutExperience[] = Array.isArray(raw.experience)
    ? raw.experience.map((it: Record<string, unknown>) => ({
        year: typeof it?.year === 'string' ? it.year : '',
        title: typeof it?.title === 'string' ? it.title : '',
        desc: typeof it?.desc === 'string' ? it.desc : '',
      }))
    : []

  const section: AboutSection[] = Array.isArray(raw.section)
    ? raw.section.map((s: Record<string, unknown>) => {
        const title = typeof s?.title === 'string' ? s.title : ''
        const items: AboutSectionItem[] = Array.isArray(s?.items)
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

  const contacts = (Array.isArray(raw.contacts) ? raw.contacts : []) as AboutContact[]

  return {
    introduction: intro,
    experience,
    section,
    contacts,
  }
}

function generateAboutJson(locale = 'zh-CN') {
  generateJsonFromYaml({
    locale,
    baseName: 'about',
    fallback: {},
    normalize: (raw) => normalize(raw as Record<string, unknown>),
  })
}

if (isDirectRun(import.meta)) {
  runCliScript('about.json', () => {
    generateAboutJson('zh-CN')
    generateAboutJson('en-US')
  })
}

export { generateAboutJson }
