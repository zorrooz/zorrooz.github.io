import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import matter from 'gray-matter';
import yaml from 'js-yaml';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentSrcDir = path.join(__dirname, '../src/content-src');
const contentOutputDir = path.join(__dirname, '../src/content');

// 确保输出目录存在
function ensureDirectoryExistence(filePath) {
  const dirname = path.dirname(filePath);
  if (fs.existsSync(dirname)) return true;
  ensureDirectoryExistence(dirname);
  fs.mkdirSync(dirname);
}

// 解析 front matter，返回 title 等元信息
function readFrontMatter(filePath) {
  try {
    const raw = fs.readFileSync(filePath, 'utf-8');
    const { data } = matter(raw);
    const fileName = path.basename(filePath, '.md');
    const title = data.title || fileName.replace(/^\d{2}-\d{2}-\d{2}--/, '');
    return {
      title,
      date: data.date || '',
      tags: Array.isArray(data.tags) ? data.tags : [],
      draft: !!data.draft,
      description: data.description || ''
    };
  } catch {
    const fileName = path.basename(filePath, '.md');
    return { title: fileName.replace(/^\d{2}-\d{2}-\d{2}--/, ''), date: '', tags: [], draft: false, description: '' };
  }
}

// 构建目录树节点（递归）
function buildNotesTree(absDir, relDirFromSrc) {
  const items = fs.readdirSync(absDir, { withFileTypes: true });
  // 只展示目录节点的 children，文件以 type: 'file' 的节点放入 children
  const node = { name: path.basename(absDir), type: 'directory', children: [] };

  for (const item of items) {
    const absPath = path.join(absDir, item.name);
    const relPathFromSrc = path.posix.join(relDirFromSrc, item.name); // 使用 posix 保持 / 分隔
    if (item.isDirectory()) {
      const childDirNode = buildNotesTree(absPath, relPathFromSrc);
      node.children.push(childDirNode);
    } else if (item.isFile() && item.name.endsWith('.md')) {
      const meta = readFrontMatter(absPath);
      // 路径以 "src/" 前缀存储，符合您的要求
      const storedPath = path.posix.join('src', 'content-src', relPathFromSrc);
      node.children.push({
        type: 'file',
        title: meta.title,
        path: storedPath,
        date: meta.date,
        tags: meta.tags,
        draft: meta.draft,
        description: meta.description
      });
    }
  }

  // 排序：目录在前，随后按名称/标题排序
  node.children.sort((a, b) => {
    const typeA = a.type === 'directory' ? 0 : 1;
    const typeB = b.type === 'directory' ? 0 : 1;
    if (typeA !== typeB) return typeA - typeB;
    const nameA = a.type === 'directory' ? a.name : a.title;
    const nameB = b.type === 'directory' ? b.name : b.title;
    return nameA.localeCompare(nameB);
  });

  return node;
}

/**
 * 生成 src/content/notes/notes.json
 * 来源：src/content-src/notes（递归）
 * 结构：{ notes: [ ...目录节点 ] }
 */
function generateNotesJson() {
  const notesSrcRoot = path.join(contentSrcDir, 'notes');
  if (!fs.existsSync(notesSrcRoot)) {
    console.log(`Source notes directory does not exist: ${notesSrcRoot}. Nothing to do.`);
    return;
  }

  // 顶层 children 作为 notes 数组
  const rootTree = buildNotesTree(notesSrcRoot, 'notes');
  const notesArray = rootTree.children;

  const outputPath = path.join(contentOutputDir, 'notes', 'notes.json');
  ensureDirectoryExistence(outputPath);
  fs.writeFileSync(outputPath, JSON.stringify({ notes: notesArray }, null, 2));
  console.log(`Successfully generated: ${outputPath}`);
}

/**
 * 将 src/content-src/notes 下的所有 Markdown 同步到 src/content/notes，
 * 并统一 Front Matter 字段结构为：
 * title, date, author, tags, draft, description
 */
function syncNotesMarkdown() {
  const notesSrcRoot = path.join(contentSrcDir, 'notes');
  const notesOutRoot = path.join(contentOutputDir, 'notes');
  if (!fs.existsSync(notesSrcRoot)) {
    console.log(`Source notes directory does not exist: ${notesSrcRoot}. Nothing to do.`);
    return;
  }

  const walk = (absDir, relDirFromNotes) => {
    const entries = fs.readdirSync(absDir, { withFileTypes: true });
    for (const ent of entries) {
      const absPath = path.join(absDir, ent.name);
      const relPath = path.join(relDirFromNotes, ent.name);
      if (ent.isDirectory()) {
        walk(absPath, relPath);
      } else if (ent.isFile() && ent.name.endsWith('.md')) {
        // 读取源 md
        const raw = fs.readFileSync(absPath, 'utf-8');
        const parsed = matter(raw);
        const data = parsed.data || {};
        const body = parsed.content || '';

        // 标准化 Front Matter
        const fileName = path.basename(absPath, '.md');
        const normalized = {
          title: data.title || fileName.replace(/^\d{2}-\d{2}-\d{2}--/, ''),
          date: data.date || '',
          author: data.author || '',
          tags: Array.isArray(data.tags) ? data.tags : [],
          draft: !!data.draft,
          description: data.description || ''
        };

        // 写入到 content/notes
        const outPath = path.join(notesOutRoot, relPath);
        ensureDirectoryExistence(outPath);
        // 生成与 r-demo 一致的多行 YAML Front Matter（双引号 + 内联数组）
        const quote = (s) => `"${String(s).replace(/"/g, '\\"')}"`;
        const arrInline = (arr) => `[${arr.map(v => quote(v)).join(', ')}]`;
        const fmLines = [
          `title: ${quote(normalized.title || '')}`,
          `date: ${quote(normalized.date || '')}`,
          `author: ${quote(normalized.author || '')}`,
          `tags: ${arrInline(Array.isArray(normalized.tags) ? normalized.tags : [])}`,
          `draft: ${normalized.draft ? 'true' : 'false'}`,
          `description: ${quote(normalized.description || '')}`
        ];
        const fm = ['---', ...fmLines, '---', '', (body || '').trim(), ''].join('\
');
        fs.writeFileSync(outPath, fm, 'utf-8');
        console.log(`Synced MD: ${outPath}`);
      }
    }
  };

  walk(notesSrcRoot, '');
}

// 主执行函数
function main() {
  console.log('Starting content generation script...');
  generateNotesJson();
  syncNotesMarkdown();
  console.log('Content generation complete.');
}

main();