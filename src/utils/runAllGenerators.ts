/**
 * Unified generator runner.
 * Order matters:
 * 1) notes, projects, topics
 * 2) categories (depends on notes/projects/topics)
 * 3) posts (depends on notes, optionally categories)
 * 4) tags (depends on posts)
 *
 * About/resources are YAML-only generators, no dependencies; run first.
 *
 * Exported as `runAllGenerators` for reuse (dev plugin, tests).
 * Auto-runs only when executed directly via `node runAllGenerators.ts`.
 */

import path from 'path'
import { fileURLToPath } from 'url'

import { generateAboutJson } from './generators/generateAbout.ts'
import { generateResourcesJson } from './generators/generateResources.ts'
import { generateNotesJson } from './generators/generateNotes.ts'
import { generateProjectsJson } from './generators/generateProjects.ts'
import { generateTopicsJson } from './generators/generateTopics.ts'
import { generateCategoriesJson } from './generators/generateCategories.ts'
import { generatePostsJson } from './generators/generatePosts.ts'
import { generateTagsJson, checkTagsConsistency } from './generators/generateTags.ts'
import { generateSitemap } from './generators/generateSitemap.ts'
import { generateHtml } from './generators/generateHtml.ts'
import { generateSearchIndex } from './generators/generateSearchIndex.ts'

const isDirectRun =
  process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)

const contentSrcDir = path.join(path.dirname(fileURLToPath(import.meta.url)), '../content-src')

async function runStep(name: string, fn: () => Promise<void> | void) {
  try {
    await fn()
  } catch (err) {
    console.error(`[Generator][${name}] failed:`, err)
    throw err
  }
}

export async function runAllGenerators() {
  console.log('== Generators: start ==')

  // 0. standalone YAML sources - Chinese + English
  await runStep('about', () => {
    generateAboutJson('zh-CN')
    generateAboutJson('en-US')
  })
  await runStep('resources', () => {
    generateResourcesJson('zh-CN')
    generateResourcesJson('en-US')
  })

  // 1. basic indexes - Chinese
  await runStep('notes', () => generateNotesJson('zh-CN'))
  await runStep('projects', () => generateProjectsJson('zh-CN'))
  await runStep('topics', () => generateTopicsJson('zh-CN'))

  // 2. categories depends on notes/projects/topics - Chinese
  await runStep('categories', () => generateCategoriesJson('zh-CN'))

  // 3. posts depends on notes - Chinese
  await runStep('posts', () => generatePostsJson('zh-CN'))

  // 4. tags depends on posts - Chinese
  await runStep('tags', () => generateTagsJson('zh-CN'))

  // 4.5 incremental translation (notes articles + about/categories/resources yaml).
  //      English yaml is generated from the Chinese yaml by the translator —
  //      hand-written -en.yaml files are not maintained. Skips up-to-date -en files;
  //      failure warns without blocking the build.
  await runStep('translate', async () => {
    if (process.env.GBLOG_NO_TRANSLATE === '1') {
      console.log('== Translators: skipped (GBLOG_NO_TRANSLATE=1) ==')
      return
    }
    try {
      const { default: manager } = await import('./translator/translator.ts')
      await manager.translate(contentSrcDir, { skipExisting: true })
    } catch (err) {
      console.warn('[Warn] incremental translation failed (build continues):', err)
    }
  })

  // 4.6 tagMerger mapping — 增量补齐 zh→en 标签映射（缓存命中 0 token）
  await runStep('tag-mapping', async () => {
    if (process.env.GBLOG_NO_TRANSLATE === '1') return
    try {
      const { ensureTagTranslation } = await import('./tagMerger/tagMerger.ts')
      await ensureTagTranslation()
    } catch (err) {
      console.warn('[Warn] tag mapping sync failed (en 标签将保持中文):', err)
    }
  })

  // 4.7 tag 一致性自动解决（tagMerger）：以 zh 文件为基准重写 -en 文件 tags
  await runStep('tag-consistency', async () => {
    if (process.env.GBLOG_NO_TRANSLATE === '1') return
    try {
      const { fixTagConsistency } = await import('./tagMerger/tagMerger.ts')
      fixTagConsistency()
    } catch (err) {
      console.warn('[Warn] tag consistency fix failed:', err)
    }
  })

  // 5. basic indexes - English
  await runStep('notes-en', () => generateNotesJson('en-US'))
  await runStep('projects-en', () => generateProjectsJson('en-US'))
  await runStep('topics-en', () => generateTopicsJson('en-US'))

  // 6. categories depends on notes/projects/topics - English
  await runStep('categories-en', () => generateCategoriesJson('en-US'))

  // 7. posts depends on notes - English
  await runStep('posts-en', () => generatePostsJson('en-US'))

  // 8. tags depends on posts - English
  await runStep('tags-en', () => generateTagsJson('en-US'))

  // 8.5 中英标签一致性校验（数量与名称必须一致）
  await runStep('tags-consistency', () => {
    const contentDir = path.join(contentSrcDir, '../content')
    checkTagsConsistency(contentDir)
  })

  // 9. article HTML (build-time markdown rendering) - Chinese + English
  await runStep('html', () => generateHtml('zh-CN'))
  await runStep('html-en', () => generateHtml('en-US'))

  // 10. sitemap + robots.txt (depends on zh categories.json)
  await runStep('sitemap', () => generateSitemap())

  // 11. search index (depends on html + categories) - Chinese + English
  await runStep('search-index', () => {
    generateSearchIndex('zh-CN')
    generateSearchIndex('en-US')
  })

  console.log('== Generators: done ==')
}

if (isDirectRun) {
  runAllGenerators().catch((err) => {
    console.error('Generators main failed:', err)
    process.exitCode = 1
  })
}
