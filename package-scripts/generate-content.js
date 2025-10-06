import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import matter from 'gray-matter';
import { marked } from 'marked';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const SRC_DIR = path.resolve(__dirname, '../src/content-src');
const OUT_DIR = path.resolve(__dirname, '../src/content');

function ensureDirSync(dirPath) {
  if (!fs.existsSync(dirPath)) fs.mkdirSync(dirPath, { recursive: true });
}

function readFileSyncSafe(filePath) {
  try {
    return fs.readFileSync(filePath, 'utf-8');
  } catch {
    return null;
  }
}

function writeJsonSync(filePath, data) {
  ensureDirSync(path.dirname(filePath));
  fs.writeFileSync(filePath, JSON.stringify(data, null, 2), 'utf-8');
}

function walkDir(dir, exts = ['.md', '.yaml', '.yml']) {
  const results = [];
  const entries = fs.readdirSync(dir, { withFileTypes: true });
  for (const entry of entries) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      results.push(...walkDir(full, exts));
    } else {
      const ext = path.extname(entry.name).toLowerCase();
      if (exts.includes(ext)) results.push(full);
    }
  }
  return results;
}

function toSlug(name) {
  return name
    .toLowerCase()
    .replace(/[^a-z0-9\u4e00-\u9fa5\-]+/gi, '-')
    .replace(/-+/g, '-')
    .replace(/^-|-$|\.md$|\.yaml$|\.yml$/g, '');
}

function firstParagraph(md) {
  const lines = md.split(/\r?\n/);
  let buf = [];
  for (const line of lines) {
    if (line.trim() === '') {
      if (buf.length) break;
      else continue;
    }
    buf.push(line.trim());
  }
  const text = buf.join(' ');
  return text.replace(/[`*_#<>]/g, '').trim();
}

function inferTypeAndCategory(relPath) {
  // relPath like notes/Bioinformatics/alignment/file.md
  const parts = relPath.split(/[\\/]/).filter(Boolean);
  const type = parts[0] || '';
  const category = parts.slice(1, Math.max(1, parts.length - 1));
  return { type, category };
}

function buildPreview(data, content) {
  const base = (data.description || data.preview || firstParagraph(content) || '').trim();
  return base.length > 200 ? base.slice(0, 197) + '...' : base;
}

function mdToHtml(md) {
  return marked.parse(md);
}

async function parseYamlContent(text) {
  try {
    const yamlMod = await import('yaml');
    return yamlMod.parse(text);
  } catch (e) {
    console.warn('YAML 解析失败，跳过：', e?.message || e);
    return null;
  }
}

function normalizeItem({ id, title, date, tags, category, preview, outRelPath }) {
  return {
    id,
    title: title || '',
    date: date || '',
    category: Array.isArray(category) ? category : category ? [category] : [],
    tags: Array.isArray(tags) ? tags : tags ? [tags] : [],
    preview: preview || '',
    path: outRelPath.replace(/\\/g, '/')
  };
}

function sortByDateDesc(a, b) {
  const da = new Date(a.date || 0).getTime();
  const db = new Date(b.date || 0).getTime();
  return db - da || (b.title || '').localeCompare(a.title || '');
}

async function main() {
  ensureDirSync(OUT_DIR);

  const files = walkDir(SRC_DIR);

  const posts = [];     // notes
  const projects = [];  // projects
  const topics = [];    // topics
  const tagCount = new Map();

  let idCounter = 1;

  for (const absPath of files) {
    const relFromSrc = path.relative(SRC_DIR, absPath);
    const ext = path.extname(absPath).toLowerCase();

    if (ext === '.md') {
      const raw = readFileSyncSafe(absPath);
      if (!raw) continue;

      const { data, content } = matter(raw);
      const { type, category } = inferTypeAndCategory(relFromSrc);
      const title = data.title || (content.match(/^#\s+(.+)$/m)?.[1] || toSlug(path.basename(absPath)));
      const date = data.date || '';
      const tags = data.tags || [];
      const preview = buildPreview(data, content);

      const html = mdToHtml(content);

      // output JSON path mirror, change extension to .json
      const outRelPath = relFromSrc.replace(/\.md$/i, '.json');
      const outAbsPath = path.join(OUT_DIR, outRelPath);
      ensureDirSync(path.dirname(outAbsPath));

      const singleJson = {
        meta: {
          title,
          date,
          author: data.author || '',
          tags,
          category,
          draft: !!data.draft,
          slug: toSlug(path.basename(absPath))
        },
        contentHtml: html
      };
      writeJsonSync(outAbsPath, singleJson);

      // aggregate item
      const item = normalizeItem({
        id: idCounter++,
        title,
        date,
        tags,
        category,
        preview,
        outRelPath
      });

      // collect into buckets
      if (type === 'notes') posts.push(item);
      else if (type === 'projects') projects.push(item);
      else if (type === 'topics') topics.push(item);

      // tag counts
      for (const t of item.tags) {
        const k = String(t);
        tagCount.set(k, (tagCount.get(k) || 0) + 1);
      }
    } else if (ext === '.yaml' || ext === '.yml') {
      const raw = readFileSyncSafe(absPath);
      if (!raw) continue;
      const parsed = await parseYamlContent(raw);
      if (!parsed) continue;

      const outRelPath = relFromSrc.replace(/\.(yaml|yml)$/i, '.json');
      const outAbsPath = path.join(OUT_DIR, outRelPath);
      writeJsonSync(outAbsPath, parsed);
    }
  }

  // sort aggregates by date desc
  posts.sort(sortByDateDesc);
  projects.sort(sortByDateDesc);
  topics.sort(sortByDateDesc);

  // write aggregate JSON
  writeJsonSync(path.join(OUT_DIR, 'posts.json'), posts);
  writeJsonSync(path.join(OUT_DIR, 'projects.json'), projects);
  writeJsonSync(path.join(OUT_DIR, 'topics.json'), topics);

  // tags.json
  const tagsArr = Array.from(tagCount.entries()).map(([name, count]) => ({ name, count }));
  writeJsonSync(path.join(OUT_DIR, 'tags.json'), tagsArr);

  console.log('内容生成完成：');
  console.log(`- 单篇 JSON 已输出至: ${OUT_DIR}/*/*.json`);
  console.log(`- 聚合: posts.json, projects.json, topics.json, tags.json`);
}

main().catch(err => {
  console.error('生成内容失败:', err);
  process.exit(1);
});