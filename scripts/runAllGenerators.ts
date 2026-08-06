/**
 * 统一生成器入口。步骤顺序有依赖关系：notes/projects/topics → categories → posts → tags，
 * 中英双 locale 各跑一轮；about/resources 为纯 YAML 独立生成，最先执行。
 * 供 dev 插件与 CLI 复用；仅直接执行时才自动运行。
 */

import { contentSrcDir, contentDir } from './dataConfig.ts'
import { isDirectRun } from './generators/core/index.ts'
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

  // 0. 独立 YAML 源 — 中英各一次
  await runStep('about', () => {
    generateAboutJson('zh-CN')
    generateAboutJson('en-US')
  })
  await runStep('resources', () => {
    generateResourcesJson('zh-CN')
    generateResourcesJson('en-US')
  })

  // 1. 基础索引 — 中文
  await runStep('notes', () => generateNotesJson('zh-CN'))
  await runStep('projects', () => generateProjectsJson('zh-CN'))
  await runStep('topics', () => generateTopicsJson('zh-CN'))

  // 2. categories 依赖 notes/projects/topics
  await runStep('categories', () => generateCategoriesJson('zh-CN'))

  // 3. posts 依赖 notes
  await runStep('posts', () => generatePostsJson('zh-CN'))

  // 4. tags 依赖 posts
  await runStep('tags', () => generateTagsJson('zh-CN'))

  // 4.5 增量翻译（notes 文章 + about/categories/resources yaml）。
  // 英文 yaml 由中文 yaml 经翻译生成，不维护手写 -en.yaml；跳过未变更的 -en 文件，
  // 失败仅告警不阻断构建。
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

  // 4.6 tagMerger 映射 — 增量补齐 zh→en 标签映射（缓存命中 0 token）
  await runStep('tag-mapping', async () => {
    if (process.env.GBLOG_NO_TRANSLATE === '1') return
    try {
      const { ensureTagTranslation } = await import('./tagMerger/tagMerger.ts')
      await ensureTagTranslation()
    } catch (err) {
      console.warn('[Warn] tag mapping sync failed (en 标签将保持中文):', err)
    }
  })

  // 4.7 tag 一致性自动解决：以 zh 文件为基准重写 -en 文件 tags
  await runStep('tag-consistency', async () => {
    if (process.env.GBLOG_NO_TRANSLATE === '1') return
    try {
      const { fixTagConsistency } = await import('./tagMerger/tagMerger.ts')
      fixTagConsistency()
    } catch (err) {
      console.warn('[Warn] tag consistency fix failed:', err)
    }
  })

  // 5. 基础索引 — 英文
  await runStep('notes-en', () => generateNotesJson('en-US'))
  await runStep('projects-en', () => generateProjectsJson('en-US'))
  await runStep('topics-en', () => generateTopicsJson('en-US'))

  // 6. categories 依赖 notes/projects/topics
  await runStep('categories-en', () => generateCategoriesJson('en-US'))

  // 7. posts 依赖 notes
  await runStep('posts-en', () => generatePostsJson('en-US'))

  // 8. tags 依赖 posts
  await runStep('tags-en', () => generateTagsJson('en-US'))

  // 8.5 中英标签一致性校验（数量与名称必须一致）
  await runStep('tags-consistency', () => {
    checkTagsConsistency(contentDir)
  })

  // 9. 文章 HTML（构建期 markdown 渲染）— 中英各一次
  await runStep('html', () => generateHtml('zh-CN'))
  await runStep('html-en', () => generateHtml('en-US'))

  // 10. sitemap + robots.txt（依赖 zh categories.json）
  await runStep('sitemap', () => generateSitemap())

  // 11. 搜索索引（依赖 html + categories）— 中英各一次
  await runStep('search-index', () => {
    generateSearchIndex('zh-CN')
    generateSearchIndex('en-US')
  })

  console.log('== Generators: done ==')
}

if (isDirectRun(import.meta)) {
  runAllGenerators().catch((err) => {
    console.error('Generators main failed:', err)
    process.exitCode = 1
  })
}
