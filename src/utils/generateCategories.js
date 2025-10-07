import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import yaml from 'js-yaml';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const contentDir = path.join(__dirname, '../content');
const contentSrcDir = path.join(__dirname, '../content-src');

const notesJsonPath = path.join(contentDir, 'notes.json');
const postsJsonPath = path.join(contentDir, 'posts.json');
const projectsJsonPath = path.join(contentDir, 'projects.json');
const topicsJsonPath = path.join(contentDir, 'topics.json');
const categoriesYamlPath = path.join(contentSrcDir, 'categories.yaml');
const categoriesJsonPath = path.join(contentDir, 'categories.json');

function ensureDirectoryExistence(filePath) {
  const dir = path.dirname(filePath);
  if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
}

function safeArray(x) {
  return Array.isArray(x) ? x : [];
}

function readJsonArray(p, fallback = []) {
  try {
    if (!fs.existsSync(p)) return fallback;
    const raw = fs.readFileSync(p, 'utf-8').trim();
    const arr = JSON.parse(raw);
    return Array.isArray(arr) ? arr : fallback;
  } catch {
    return fallback;
  }
}

function readYaml(p, fallback = {}) {
  try {
    if (!fs.existsSync(p)) return fallback;
    const raw = fs.readFileSync(p, 'utf-8');
    const obj = yaml.load(raw);
    return obj || fallback;
  } catch {
    return fallback;
  }
}

function normalizeNoteYaml(it) {
  const categories = it?.categories && typeof it.categories === 'object' ? it.categories : {};
  return {
    name: typeof it?.name === 'string' ? it.name : '',
    title: typeof it?.title === 'string' ? it.title : '',
    desc: typeof it?.desc === 'string' ? it.desc : '',
    date: typeof it?.date === 'string' ? it.date : '',
    categories
  };
}

function normalizeProjectOrTopicYaml(it) {
  const categories = it?.categories && typeof it.categories === 'object' ? it.categories : {};
  return {
    name: typeof it?.name === 'string' ? it.name : '',
    title: typeof it?.title === 'string' ? it.title : '',
    desc: typeof it?.desc === 'string' ? it.desc : '',
    date: typeof it?.date === 'string' ? it.date : '',
    tags: Array.isArray(it?.tags) ? it.tags.filter(t => typeof t === 'string') : [],
    github: typeof it?.github === 'string' ? it.github : '',
    doi: typeof it?.doi === 'string' ? it.doi : '',
    categories
  };
}

function deriveGroupKey(post) {
  // posts.json: category [Group, Subgroup]
  const cat = safeArray(post?.category);
  if (cat.length >= 1) return cat[0];
  return 'Other';
}

function aggregateNotesByYaml(noteDefs, postsArr, notesArr) {
  // title -> notes.relativePath, 供拼装 articleUrl
  const noteByTitle = new Map();
  for (const n of notesArr) {
    if (typeof n?.title === 'string' && typeof n?.relativePath === 'string') {
      noteByTitle.set(n.title, n.relativePath);
    }
  }

  const items = [];

  for (const defRaw of noteDefs) {
    const def = normalizeNoteYaml(defRaw);
    if (!def.name) continue;

    // 使用 name 作为 posts 的一级分类匹配键
    const groupKey = def.name;
    const arr = postsArr.filter(p => deriveGroupKey(p) === groupKey);

    // 统计
    const postsCount = arr.length;
    const words = arr.reduce((sum, x) => sum + (Number(x?.wordCount) || 0), 0);
    const latest = arr.map(x => Date.parse(x?.date || '') || 0).reduce((a, b) => Math.max(a, b), 0);
    const latestDate = def.date || (latest ? new Date(latest).toISOString().slice(0, 10) : '');

    // tags 去重
    const tagSet = new Set();
    for (const x of arr) {
      safeArray(x?.tags).forEach(t => {
        if (typeof t === 'string' && t.trim()) tagSet.add(t.trim());
      });
    }
    const tags = Array.from(tagSet);

    // 文章列表
    const articles = arr.map(x => {
      const rel = noteByTitle.get(x.title);
      const articleUrl = rel ? `/article/notes/${rel}` : '';
      return {
        title: x.title || '',
        articleUrl
      };
    });

    // root：第一篇（按 relativePath 字典序）文章
    const rels = notesArr
      .filter(n => Array.isArray(n?.relativePath ? [n.relativePath] : []) && arr.some(p => p.title === n.title))
      .map(n => n.relativePath)
      .filter(Boolean)
      .sort((a, b) => a.localeCompare(b));

    const root = rels.length ? `/article/notes/${rels[0]}` : (articles.find(a => a.articleUrl)?.articleUrl || '');

    items.push({
      name: def.name,
      title: def.title || def.name,
      desc: def.desc,
      date: latestDate,
      tags,
      postsCount,
      words,
      articles,
      root,
      categories: def.categories
    });
  }

  // 稳定排序：按名称
  items.sort((a, b) => a.name.localeCompare(b.name, 'zh-Hans-CN'));
  return items;
}

function normalizeYamlRoot(yobj) {
  return {
    notes: Array.isArray(yobj?.notes) ? yobj.notes : [],
    projects: Array.isArray(yobj?.projects) ? yobj.projects : [],
    topics: Array.isArray(yobj?.topics) ? yobj.topics : []
  };
}

function buildProjectOrTopicItems(defs, contentArr, kind) {
  const items = [];

  for (const defRaw of defs) {
    const def = normalizeProjectOrTopicYaml(defRaw);
    if (!def.name) continue;

    // 以 name 的小写去空格形式匹配；再加入“仅字母数字”的归一化键；若 url 是 /article/ 则解析首段加入候选
    const nameKey = def.name.toLowerCase().replace(/\s+/g, '');
    const nameAlpha = def.name.toLowerCase().replace(/[^a-z0-9]+/g, ''); // 例："Demo 课题" -> "demo"
    let urlFirstSeg = '';
    if (typeof def.url === 'string' && def.url.startsWith('/article/')) {
      const segs = def.url.replace(/^\/+/, '').split('/');
      if (segs.length >= 2) urlFirstSeg = segs[1].toLowerCase();
    }
    const candidates = Array.from(new Set([nameKey, nameAlpha, urlFirstSeg].filter(Boolean)));

    const arr = contentArr.filter(x => {
      if (typeof x?.relativePath !== 'string') return false;
      const firstRaw = x.relativePath.split('/')[0];
      const firstA = firstRaw.toLowerCase().replace(/\s+/g, '');
      const firstAlpha = firstRaw.toLowerCase().replace(/[^a-z0-9]+/g, '');
      return candidates.includes(firstA) || candidates.includes(firstAlpha);
    });

    const postsCount = arr.length;
    const words = arr.reduce((s, x) => s + (Number(x?.wordCount) || 0), 0);

    const articles = arr.map(n => ({
      title: typeof n.title === 'string' ? n.title : '',
      articleUrl: `/article/${kind === 'project' ? 'projects' : 'topics'}/${n.relativePath}`
    }));

    // root：第一篇（相对路径字典序）；若无匹配且有 def.url，则用 def.url 兜底
    const rels = arr.map(n => n.relativePath).filter(Boolean).sort((a, b) => a.localeCompare(b));
    const root = rels.length
      ? `/article/${kind === 'project' ? 'projects' : 'topics'}/${rels[0]}`
      : (def.url && def.url.startsWith('/article/') ? def.url : '');

    const out = {
      name: def.name,
      title: def.title || def.name,
      desc: def.desc,
      date: def.date || '',
      tags: def.tags || [],
      postsCount,
      words,
      articles,
      root,
      categories: def.categories
    };

    // 附加外链字段
    if (kind === 'project' && def.github) out.github = def.github;
    if (kind === 'topic' && def.doi) out.doi = def.doi;

    items.push(out);
  }

  // 稳定排序：按名称
  items.sort((a, b) => a.name.localeCompare(b.name, 'zh-Hans-CN'));
  return items;
}

function mergeSections(notesItems, projectItems, topicItems) {
  return [
    { title: '笔记', items: notesItems },
    { title: '项目', items: projectItems },
    { title: '课题', items: topicItems }
  ];
}

function generateCategoriesJson() {
  // 输入
  const notesArr = readJsonArray(notesJsonPath, []);
  const postsArr = readJsonArray(postsJsonPath, []);
  const projectsArr = readJsonArray(projectsJsonPath, []);
  const topicsArr = readJsonArray(topicsJsonPath, []);
  const yamlObj = readYaml(categoriesYamlPath, {});
  const { notes: noteDefs, projects: projectDefs, topics: topicDefs } = normalizeYamlRoot(yamlObj);

  // 构建各区块
  const notesItems = aggregateNotesByYaml(noteDefs, postsArr, notesArr);
  const projectItems = buildProjectOrTopicItems(projectDefs, projectsArr, 'project');
  const topicItems = buildProjectOrTopicItems(topicDefs, topicsArr, 'topic');

  // 合并区块（直接输出数组，移除 categoryList 包裹）
  const sections = mergeSections(notesItems, projectItems, topicItems);

  try {
    ensureDirectoryExistence(categoriesJsonPath);
    fs.writeFileSync(categoriesJsonPath, JSON.stringify(sections, null, 2), 'utf-8');
    console.log(`Successfully generated: ${categoriesJsonPath}`);
  } catch (e) {
    console.error('Failed to write categories.json:', e?.message || e);
    process.exitCode = 1;
  }
}

function main() {
  console.log('Starting categories.json generation script...');
  generateCategoriesJson();
  console.log('categories.json generation complete.');
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main();
}

export { generateCategoriesJson };