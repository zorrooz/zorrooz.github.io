import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentDir = path.join(__dirname, '../content');
const postsJsonPath = path.join(contentDir, 'posts.json');
const tagsJsonPath = path.join(contentDir, 'tags.json');

function ensureDirectoryExistence(filePath) {
  const dir = path.dirname(filePath);
  if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
}

function generateTagsFromPosts(posts) {
  const map = new Map();
  for (const p of posts) {
    const tags = Array.isArray(p?.tags) ? p.tags : [];
    for (const raw of tags) {
      const name = typeof raw === 'string' ? raw.trim() : '';
      if (!name) continue;
      map.set(name, (map.get(name) || 0) + 1);
    }
  }
  // 输出为 { name, count }，按名称排序，便于组件渲染稳定
  return Array.from(map.entries())
    .map(([name, count]) => ({ name, count }))
    .sort((a, b) => a.name.localeCompare(b.name, 'zh-Hans-CN'));
}

function generateTagsJson() {
  if (!fs.existsSync(postsJsonPath)) {
    console.error(`posts.json not found at ${postsJsonPath}. Please generate posts.json first.`);
    process.exitCode = 1;
    return;
  }

  let posts = [];
  try {
    const raw = fs.readFileSync(postsJsonPath, 'utf-8').trim();
    posts = JSON.parse(raw);
    if (!Array.isArray(posts)) {
      console.error('posts.json is not an array. Abort.');
      process.exitCode = 1;
      return;
    }
  } catch (e) {
    console.error('Failed to read/parse posts.json:', e?.message || e);
    process.exitCode = 1;
    return;
  }

  const tags = generateTagsFromPosts(posts);

  try {
    ensureDirectoryExistence(tagsJsonPath);
    fs.writeFileSync(tagsJsonPath, JSON.stringify(tags, null, 2), 'utf-8');
    console.log(`Successfully generated: ${tagsJsonPath} (${tags.length} tags)`);
  } catch (e) {
    console.error('Failed to write tags.json:', e?.message || e);
    process.exitCode = 1;
  }
}

function main() {
  console.log('Starting tags.json generation script...');
  generateTagsJson();
  console.log('tags.json generation complete.');
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main();
}

export { generateTagsJson };