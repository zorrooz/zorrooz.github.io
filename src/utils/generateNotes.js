import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import yaml from 'js-yaml';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentSrcDir = path.join(__dirname, '../content-src');
const notesSrcDir = path.join(contentSrcDir, 'notes');
const contentOutputDir = path.join(__dirname, '../content');

// 根据语言获取对应的文件路径
function getFilePaths(locale = 'zh-CN') {
  const suffix = locale === 'zh-CN' ? '' : '-en';
  return {
    outputPath: path.join(contentOutputDir, `notes${suffix}.json`)
  };
}

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

function normalizeTags(tags) {
  if (Array.isArray(tags)) {
    return tags
      .map(t => (typeof t === 'string' ? t.trim() : ''))
      .filter(Boolean);
  }
  if (typeof tags === 'string') {
    // support comma/space separated
    return tags
      .split(/[,，]/g)
      .map(s => s.trim())
      .filter(Boolean);
  }
  return [];
}

function parseFrontMatterAndBody(raw) {
  // Support YAML frontmatter delimited by --- at start of file
  if (raw.startsWith('---')) {
    const lines = raw.split(/\r?\n/);
    let endIdx = -1;
    for (let i = 1; i < lines.length; i++) {
      if (lines[i].trim() === '---') {
        endIdx = i;
        break;
      }
    }
    if (endIdx !== -1) {
      const fmText = lines.slice(1, endIdx).join('\n');
      let fm = {};
      try {
        fm = yaml.load(fmText) || {};
      } catch (e) {
        console.warn('Warn: failed to parse frontmatter YAML. Using empty object.', e?.message || e);
      }
      const body = lines.slice(endIdx + 1).join('\n');
      return { frontmatter: fm || {}, body };
    }
  }
  return { frontmatter: {}, body: raw };
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

function buildNoteItem(mdPath) {
  const raw = fs.readFileSync(mdPath, 'utf-8');
  const { frontmatter, body } = parseFrontMatterAndBody(raw);

  const title = typeof frontmatter?.title === 'string' ? frontmatter.title : '';
  const date = typeof frontmatter?.date === 'string' ? frontmatter.date : '';
  const description = typeof frontmatter?.description === 'string' ? frontmatter.description : '';
  const draft = typeof frontmatter?.draft === 'boolean' ? frontmatter.draft : false;
  const tags = normalizeTags(frontmatter?.tags);

  const plain = markdownToPlain(body);
  const wordCount = countWordsSmart(plain);
  const tagCount = tags.length;

  const relativePath = toPosixRelativeNoExt(mdPath, notesSrcDir);

  return {
    title,
    date,
    tags,
    draft,
    description,
    relativePath,
    wordCount,
    tagCount
  };
}

function generateNotesJson(locale = 'zh-CN') {
  const paths = getFilePaths(locale);
  
  if (!fs.existsSync(notesSrcDir)) {
    console.warn(`Warn: notes source directory not found at ${notesSrcDir}`);
  }
  
  // 根据语言过滤对应的markdown文件
  const mdFiles = walk(notesSrcDir, p => {
    if (locale === 'zh-CN') {
      // 中文版本：只包含非英文文件（不以-en.md结尾）
      return /\.md$/i.test(p) && !p.endsWith('-en.md');
    } else {
      // 英文版本：只包含英文文件（以-en.md结尾）
      return /-en\.md$/i.test(p);
    }
  });

  const items = mdFiles.map(buildNoteItem);

  // Optional: sort by date desc if date comparable, fallback by path
  items.sort((a, b) => {
    const ad = Date.parse(a.date || '') || 0;
    const bd = Date.parse(b.date || '') || 0;
    if (ad !== bd) return bd - ad;
    return a.relativePath.localeCompare(b.relativePath);
  });

  try {
    const targetPath = paths.outputPath;
    ensureDirectoryExistence(targetPath);
    fs.writeFileSync(targetPath, JSON.stringify(items, null, 2), 'utf-8');
    console.log(`Successfully generated: ${targetPath} (${items.length} notes)`);
  } catch (e) {
    console.error(`Failed to write ${locale} notes.json:`, e);
  }
}

function main() {
  console.log('Starting notes.json generation script...');
  generateNotesJson('zh-CN');
  generateNotesJson('en-US');
  console.log('notes.json generation complete.');
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  // executed directly
  main();
}

export { generateNotesJson };