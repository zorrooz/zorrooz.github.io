import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentSrcDir = path.join(__dirname, '../content-src');
const projectsSrcDir = path.join(contentSrcDir, 'projects');
const contentOutputDir = path.join(__dirname, '../content');
const outputPath = path.join(contentOutputDir, 'projects.json');

function ensureDirectoryExistence(filePath) {
  const dirname = path.dirname(filePath);
  if (!fs.existsSync(dirname)) {
    fs.mkdirSync(dirname, { recursive: true });
  }
}

// Recursively list files under dir that match predicate
function walk(dir, predicate = () => true, acc = []) {
  if (!fs.existsSync(dir)) return acc;
  const items = fs.readdirSync(dir, { withFileTypes: true });
  for (const it of items) {
    const full = path.join(dir, it.name);
    if (it.isDirectory()) walk(full, predicate, acc);
    else if (predicate(full)) acc.push(full);
  }
  return acc;
}

// Crude markdown -> plain text for counting words
function markdownToPlain(text) {
  // remove code fences first
  let t = text.replace(/```[\s\S]*?```/g, ' ');
  // remove inline code
  t = t.replace(/`[^`]*`/g, ' ');
  // remove images ![alt](url)
  t = t.replace(/!\[[^\]]*]\([^)]+\)/g, ' ');
  // remove links [text](url)
  t = t.replace(/\[[^\]]*]\([^)]+\)/g, ' ');
  // remove html tags
  t = t.replace(/<[^>]+>/g, ' ');
  // remove headings/bold/italic/quotes/lists/tables markers
  t = t.replace(/^#{1,6}\s+/gm, ' ');
  t = t.replace(/[*_~`>#|-]{1,}/g, ' ');
  // remove reference links [id]: url
  t = t.replace(/^\s*\[[^\]]+]:\s+\S+.*$/gm, ' ');
  // collapse whitespace
  t = t.replace(/\s+/g, ' ').trim();
  return t;
}

function countWordsSmart(text) {
  // Separate CJK chars and non-CJK words
  const cjk = (text.match(/[\u4E00-\u9FFF\u3400-\u4DBF]/g) || []).length;
  const words = (text.match(/[A-Za-z0-9]+(?:'[A-Za-z0-9]+)?/g) || []).length;
  return cjk + words;
}

function toPosixRelativeNoExt(fullPath, baseDir) {
  const rel = path.relative(baseDir, fullPath);
  const noExt = rel.replace(/\.[^/.]+$/, '');
  // normalize to POSIX style
  return noExt.split(path.sep).join('/');
}

function extractTitleFromH1(raw) {
  // find first level-1 heading line: starts with "# " (but not "##")
  const lines = raw.split(/\r?\n/);
  for (const line of lines) {
    const m = line.match(/^#\s+(.*)$/);
    if (m) return m[1].trim();
  }
  return '';
}

function buildProjectItem(mdPath) {
  const raw = fs.readFileSync(mdPath, 'utf-8');
  const title = extractTitleFromH1(raw) || path.basename(mdPath).replace(/\.md$/i, '');
  const plain = markdownToPlain(raw);
  const wordCount = countWordsSmart(plain);
  const relativePath = toPosixRelativeNoExt(mdPath, projectsSrcDir);

  return { title, relativePath, wordCount };
}

function generateProjectsJson() {
  if (!fs.existsSync(projectsSrcDir)) {
    console.warn(`Warn: projects source directory not found at ${projectsSrcDir}`);
  }
  const mdFiles = walk(projectsSrcDir, p => /\.md$/i.test(p));

  const items = mdFiles.map(buildProjectItem);

  // Sort by relativePath for stable navigation
  items.sort((a, b) => a.relativePath.localeCompare(b.relativePath));

  try {
    ensureDirectoryExistence(outputPath);
    fs.writeFileSync(outputPath, JSON.stringify(items, null, 2), 'utf-8');
    console.log(`Successfully generated: ${outputPath} (${items.length} projects)`);
  } catch (e) {
    console.error('Failed to write projects.json:', e);
  }
}

function main() {
  console.log('Starting projects.json generation script...');
  generateProjectsJson();
  console.log('projects.json generation complete.');
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  // executed directly
  main();
}

export { generateProjectsJson };