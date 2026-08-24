/**
 * 统一生成器入口。步骤顺序有依赖关系：notes/projects/topics → categories → posts → tags，
 * 中英双 locale 各跑一轮；about/resources 为纯 YAML 独立生成，最先执行。
 * 供 dev 插件与 CLI 复用；仅直接执行时才自动运行。
 */

import { contentSrcDir, contentDir } from './dataConfig.ts'
import { isDirectRun } from './lib/cli.ts'
import { logError, logInfo, logSection, logWarn } from './lib/log.ts'
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
    logError(`生成器 ${name} 失败:`, err)
    throw err
  }
}

type Locale = 'zh-CN' | 'en-US'

/** 单个 locale 的基础索引步骤（顺序有依赖：notes/projects/topics → categories → posts → tags） */
const localeSteps = (locale: Locale): Array<{ name: string; run: () => void }> => {
  const suffix = locale === 'en-US' ? '-en' : ''
  return [
    { name: `notes${suffix}`, run: () => generateNotesJson(locale) },
    { name: `projects${suffix}`, run: () => generateProjectsJson(locale) },
    { name: `topics${suffix}`, run: () => generateTopicsJson(locale) },
    { name: `categories${suffix}`, run: () => generateCategoriesJson(locale) },
    { name: `posts${suffix}`, run: () => generatePostsJson(locale) },
    { name: `tags${suffix}`, run: () => generateTagsJson(locale) },
  ]
}

async function runLocaleBlock(locale: Locale) {
  for (const step of localeSteps(locale)) {
    await runStep(step.name, step.run)
  }
}

export async function runAllGenerators() {
  logSection('Generators 开始')

  /* 0. 独立 YAML 源 — 中英各一次 */
  await runStep('about', () => {
    generateAboutJson('zh-CN')
    generateAboutJson('en-US')
  })
  await runStep('resources', () => {
    generateResourcesJson('zh-CN')
    generateResourcesJson('en-US')
  })

  /* 1-4. 中文基础索引（notes/projects/topics → categories → posts → tags） */
  await runLocaleBlock('zh-CN')

  /**
   * 4.5 增量翻译（notes 文章 + about/categories/resources yaml）。
   * 英文 yaml 由中文 yaml 翻译生成，不维护手写 -en.yaml；
   * 跳过未变更文件，失败仅告警不阻断构建。
   */
  await runStep('translate', async () => {
    if (process.env.GBLOG_NO_TRANSLATE === '1') {
      logInfo('翻译步骤已跳过（GBLOG_NO_TRANSLATE=1）')
      return
    }
    try {
      const { default: manager } = await import('./tools/translator/translator.ts')
      await manager.translate(contentSrcDir, { skipExisting: true })
    } catch (err) {
      logWarn('增量翻译失败（构建继续）:', err)
    }
  })

  /* 4.6 tagMerger 映射 — 增量补齐 zh→en 标签映射（缓存命中 0 token） */
  await runStep('tag-mapping', async () => {
    if (process.env.GBLOG_NO_TRANSLATE === '1') return
    try {
      const { ensureTagTranslation } = await import('./tools/tagMerger/tagMerger.ts')
      await ensureTagTranslation()
    } catch (err) {
      logWarn('标签映射同步失败（en 标签将保持中文）:', err)
    }
  })

  /* 4.7 tag 一致性自动解决：以 zh 文件为基准重写 -en 文件 tags */
  await runStep('tag-consistency', async () => {
    if (process.env.GBLOG_NO_TRANSLATE === '1') return
    try {
      const { fixTagConsistency } = await import('./tools/tagMerger/tagMerger.ts')
      fixTagConsistency()
    } catch (err) {
      logWarn('标签一致性修复失败:', err)
    }
  })

  /* 5-8. 英文基础索引（依赖上方翻译步骤产物） */
  await runLocaleBlock('en-US')

  /* 8.5 中英标签一致性校验（数量与名称必须一致） */
  await runStep('tags-consistency', () => {
    checkTagsConsistency(contentDir)
  })

  /* 9. 文章 HTML（构建期 markdown 渲染）— 中英各一次 */
  await runStep('html', () => generateHtml('zh-CN'))
  await runStep('html-en', () => generateHtml('en-US'))

  /* 10. sitemap + robots.txt */
  await runStep('sitemap', () => generateSitemap())

  /* 11. 搜索索引 — 中英各一次 */
  await runStep('search-index', () => {
    generateSearchIndex('zh-CN')
    generateSearchIndex('en-US')
  })

  logSection('Generators 完成')
}

if (isDirectRun(import.meta)) {
  runAllGenerators().catch((err) => {
    logError('生成器主流程失败:', err)
    process.exitCode = 1
  })
}
