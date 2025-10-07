import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentDir = path.join(__dirname, '../content');
const notesJsonPath = path.join(contentDir, 'notes.json');
const postsJsonPath = path.join(contentDir, 'posts.json');
const categoriesJsonPath = path.join(contentDir, 'categories.json');

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

// 从 categories.json 构建 notes 分区的映射：groupName -> { groupTitle, subMap }
function buildNotesCategoryMap(categoriesArr) {
  const map = new Map();
  if (!Array.isArray(categoriesArr)) return map;
  // 找到“笔记”分区
  const notesSection = categoriesArr.find(s => s && s.title === '笔记' && Array.isArray(s.items));
  if (!notesSection) return map;
  for (const item of notesSection.items) {
    if (!item || typeof item !== 'object') continue;
    const name = item.name; // 目录名，如 Omics / Programming
    if (typeof name !== 'string') continue;
    const groupTitle = typeof item.title === 'string' && item.title.trim() ? item.title : name;
    const subMap = item.categories && typeof item.categories === 'object' ? item.categories : {};
    map.set(name, { groupTitle, subMap });
  }
  return map;
}

function buildPostsFromNotes(notesArr, notesCategoryMap) {
  const posts = [];
  let id = 1;

  for (const it of notesArr) {
    if (!it || typeof it !== 'object') continue;
    if (it.draft === true) continue; // 跳过草稿

    const title = typeof it.title === 'string' ? it.title : '';
    const date = typeof it.date === 'string' ? it.date : '';
    const tags = Array.isArray(it.tags) ? it.tags : [];
    const preview = typeof it.description === 'string' ? it.description : '';
    // 基于 categories.json 的映射来生成分类显示名
    const [grp, sub] = deriveCategory(it.relativePath);
    let category = [];
    if (grp) {
      const entry = notesCategoryMap && notesCategoryMap.get(grp);
      if (entry) {
        const first = entry.groupTitle || grp;
        const second = sub ? (entry.subMap?.[sub] || sub) : undefined;
        category = second ? [first, second] : [first];
      } else {
        category = sub ? [grp, sub] : [grp];
      }
    }
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

  // 读取 categories.json（若不存在则回退老逻辑）
  let categoriesArr = [];
  try {
    if (fs.existsSync(categoriesJsonPath)) {
      const rawCat = fs.readFileSync(categoriesJsonPath, 'utf-8').trim();
      categoriesArr = JSON.parse(rawCat);
    }
  } catch {
    categoriesArr = [];
  }
  const notesCategoryMap = buildNotesCategoryMap(categoriesArr);

  const posts = buildPostsFromNotes(notesArr, notesCategoryMap);

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