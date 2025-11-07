/**
 * Unified generator runner.
 * Order matters:
 * 1) notes, projects, topics
 * 2) categories (depends on notes/projects/topics)
 * 3) posts (depends on notes, optionally categories)
 * 4) tags (depends on posts)
 * 
 * About/resources execute on import (their files call main() at top-level),
 * so importing them once is enough; do NOT call their exported functions again.
 */

import './generators/generateAbout.js';       // self-executes on import
import './generators/generateResources.js';   // self-executes on import

import { generateNotesJson } from './generators/generateNotes.js';
import { generateProjectsJson } from './generators/generateProjects.js';
import { generateTopicsJson } from './generators/generateTopics.js';
import { generateCategoriesJson } from './generators/generateCategories.js';
import { generatePostsJson } from './generators/generatePosts.js';
import { generateTagsJson } from './generators/generateTags.js';

function runStep(name, fn) {
  try {
    const ret = fn();
    // handle sync or promise
    if (ret && typeof ret.then === 'function') {
      return ret.catch(err => {
        console.error(`[Generator][${name}] failed:`, err);
        process.exitCode = 1;
      });
    }
  } catch (err) {
    console.error(`[Generator][${name}] failed:`, err);
    process.exitCode = 1;
  }
}

async function main() {
  console.log('== Generators: start ==');

  // 1. basic indexes - Chinese
  await runStep('notes', () => generateNotesJson('zh-CN'));
  await runStep('projects', () => generateProjectsJson('zh-CN'));
  await runStep('topics', () => generateTopicsJson('zh-CN'));

  // 2. categories depends on notes/projects/topics - Chinese
  await runStep('categories', () => generateCategoriesJson('zh-CN'));

  // 3. posts depends on notes - Chinese
  await runStep('posts', () => generatePostsJson('zh-CN'));

  // 4. tags depends on posts - Chinese
  await runStep('tags', () => generateTagsJson('zh-CN'));

  // 5. basic indexes - English
  await runStep('notes-en', () => generateNotesJson('en-US'));
  await runStep('projects-en', () => generateProjectsJson('en-US'));
  await runStep('topics-en', () => generateTopicsJson('en-US'));

  // 6. categories depends on notes/projects/topics - English
  await runStep('categories-en', () => generateCategoriesJson('en-US'));

  // 7. posts depends on notes - English
  await runStep('posts-en', () => generatePostsJson('en-US'));

  // 8. tags depends on posts - English
  await runStep('tags-en', () => generateTagsJson('en-US'));

  console.log('== Generators: done ==');
}

main().catch(err => {
  console.error('Generators main failed:', err);
  process.exitCode = 1;
});