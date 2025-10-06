import fs from 'fs';
import path from 'path';
import matter from 'gray-matter';
import { fileURLToPath, pathToFileURL } from 'url';

/**
 * 生成 src/content/notes/notes.json
 * - 正确解析 Markdown frontmatter 元信息
 * - 保留目录树结构（children + files）
 * - 路径相对 src/content，分隔符标准化为 /
 */

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentDir = path.join(__dirname, '../content');
const notesDir = path.join(contentDir, 'notes');

/** 规范化相对路径为 posix 风格 */
function toPosixRelative(filePath) {
  return path.relative(contentDir, filePath).replace(/\\/g, '/');
}

/** 从文件解析元数据，带稳健回退 */
function getFileMetadata(filePath) {
  const raw = fs.readFileSync(filePath, 'utf-8');
  const parsed = matter(raw || '');
  const data = parsed.data || {};
  const fileName = path.basename(filePath);

  // 依据文件名生成默认标题（去掉日期前缀与扩展名）
  const defaultTitle = fileName
    .replace(/\.md$/i, '')
    .replace(/^\d{2}-\d{2}-\d{2}--/, '')
    .replace(/[-_]+/g, ' ')
    .trim();

  const meta = {
    title: data.title || defaultTitle,
    path: toPosixRelative(filePath),
  };

  // 可选字段（若存在则保留）
  if (data.date) meta.date = data.date;
  if (data.tags) meta.tags = Array.isArray(data.tags)
    ? data.tags
    : String(data.tags).split(',').map(s => s.trim()).filter(Boolean);
  if (data.summary) meta.summary = data.summary;
  if (data.slug) meta.slug = data.slug;

  return meta;
}

/** 构建目录树：目录节点含 children，文件集合放入 files */
function buildDirectoryTree(dir) {
  const items = fs.readdirSync(dir, { withFileTypes: true });

  const children = [];
  const files = [];

  for (const item of items) {
    const fullPath = path.join(dir, item.name);
    if (item.isDirectory()) {
      children.push(buildDirectoryTree(fullPath));
    } else if (item.isFile() && /\.md$/i.test(item.name)) {
      files.push(getFileMetadata(fullPath));
    }
  }

  const node = {
    name: path.basename(dir),
    type: 'directory',
  };
  if (children.length) node.children = children;
  if (files.length) node.files = files;

  return node;
}

/** 主入口：生成 notes.json */
export function generateNotes() {
  if (!fs.existsSync(notesDir)) {
    throw new Error(`Notes 目录不存在：${notesDir}`);
  }

  const treeRoot = buildDirectoryTree(notesDir);

  // 与现有使用保持一致：仅输出子节点
  const output = { notes: treeRoot.children || [] };

  const outputPath = path.join(notesDir, 'notes.json');
  fs.writeFileSync(outputPath, JSON.stringify(output, null, 2), 'utf-8');

  console.log('notes.json 已生成：', outputPath);
}

/** 直接执行文件时运行（ESM 环境下的标准判断） */
if (import.meta.url === pathToFileURL(__filename).href) {
  generateNotes();
}