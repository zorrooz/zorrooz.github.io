import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentDir = path.join(__dirname, '../content');
const notesJsonPath = path.join(contentDir, 'notes.json');
const postsJsonPath = path.join(contentDir, 'posts.json');

function ensureDirectoryExistence(filePath) {
  const dir = path.dirname(filePath);
  if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
}

function toArray(x) {
  return Array.isArray(x) ? x : (x == null ? [] : [x]);
}

function deriveCategory(relativePath) {
  if (typeof relativePath !== 'string' || !relativePath) return [];
  const parts = relativePath.split('/'); // already POSIX style
  // 期望结构: Group/Subgroup/.../slug
  if (parts.length >= 2) return [parts[0], parts[1]];
  if (parts.length === 1) return [parts[0]];
  return [];
}

function buildPostsFromNotes(notesArr) {
  const posts = [];
  let id = 1;

  for (const it of notesArr) {
    if (!it || typeof it !== 'object') continue;
    if (it.draft === true) continue; // 跳过草稿

    const title = typeof it.title === 'string' ? it.title : '';
    const date = typeof it.date === 'string' ? it.date : '';
    const tags = Array.isArray(it.tags) ? it.tags : [];
    const preview = typeof it.description === 'string' ? it.description : '';
    const category = deriveCategory(it.relativePath);
    const wordCount = Number.isFinite(it.wordCount) ? it.wordCount : 0;

    posts.push({
      id: id,
      no: id,
      title,
      date,
      category,
      tags,
      preview,
      wordCount
    });

    id += 1;
  }

  return posts;
}

function generatePostsJson() {
  if (!fs.existsSync(notesJsonPath)) {
    console.error(`notes.json not found at ${notesJsonPath}. Please run generateNotes.js first.`);
    process.exitCode = 1;
    return;
  }

  let notesArr = [];
  try {
    const raw = fs.readFileSync(notesJsonPath, 'utf-8').trim();
    // notes.json 是数组
    notesArr = JSON.parse(raw);
    if (!Array.isArray(notesArr)) {
      console.error('notes.json is not an array. Abort.');
      process.exitCode = 1;
      return;
    }
  } catch (e) {
    console.error('Failed to read/parse notes.json:', e?.message || e);
    process.exitCode = 1;
    return;
  }

  const posts = buildPostsFromNotes(notesArr);

  try {
    ensureDirectoryExistence(postsJsonPath);
    fs.writeFileSync(postsJsonPath, JSON.stringify(posts, null, 2), 'utf-8');
    console.log(`Successfully generated: ${postsJsonPath} (${posts.length} posts)`);
  } catch (e) {
    console.error('Failed to write posts.json:', e?.message || e);
    process.exitCode = 1;
  }
}

function main() {
  console.log('Starting posts.json generation script...');
  generatePostsJson();
  console.log('posts.json generation complete.');
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main();
}

export { generatePostsJson };