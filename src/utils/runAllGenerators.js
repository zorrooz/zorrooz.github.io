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
 * Auto-runs only when executed directly via `node runAllGenerators.js`.
 */

import path from 'path'
import { fileURLToPath } from 'url'

import { generateAboutJson } from './generators/generateAbout.js'
import { generateResourcesJson } from './generators/generateResources.js'
import { generateNotesJson } from './generators/generateNotes.js'
import { generateProjectsJson } from './generators/generateProjects.js'
import { generateTopicsJson } from './generators/generateTopics.js'
import { generateCategoriesJson } from './generators/generateCategories.js'
import { generatePostsJson } from './generators/generatePosts.js'
import { generateTagsJson } from './generators/generateTags.js'
import { generateHtml } from './generators/generateHtml.js'

const isDirectRun =
  process.argv[1] &&
  path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)

async function runStep(name, fn) {
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

  // 9. article HTML (build-time markdown rendering) - Chinese + English
  await runStep('html', () => generateHtml('zh-CN'))
  await runStep('html-en', () => generateHtml('en-US'))

  console.log('== Generators: done ==')
}

if (isDirectRun) {
  runAllGenerators().catch((err) => {
    console.error('Generators main failed:', err)
    process.exitCode = 1
  })
}
