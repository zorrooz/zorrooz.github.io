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

import './generateAbout.js';       // self-executes on import
import './generateResources.js';   // self-executes on import

import { generateNotesJson } from './generateNotes.js';
import { generateProjectsJson } from './generateProjects.js';
import { generateTopicsJson } from './generateTopics.js';
import { generateCategoriesJson } from './generateCategories.js';
import { generatePostsJson } from './generatePosts.js';
import { generateTagsJson } from './generateTags.js';

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

  // 1. basic indexes
  await runStep('notes', generateNotesJson);
  await runStep('projects', generateProjectsJson);
  await runStep('topics', generateTopicsJson);

  // 2. categories depends on notes/projects/topics
  await runStep('categories', generateCategoriesJson);

  // 3. posts depends on notes (and uses categories if present)
  await runStep('posts', generatePostsJson);

  // 4. tags depends on posts
  await runStep('tags', generateTagsJson);

  console.log('== Generators: done ==');
}

main().catch(err => {
  console.error('Generators main failed:', err);
  process.exitCode = 1;
});