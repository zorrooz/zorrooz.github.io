import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import yaml from 'js-yaml';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// 路径配置（与原逻辑一致）
const contentDir = path.join(__dirname, '../content');
const contentSrcDir = path.join(__dirname, '../content-src');
const notesJsonPath = path.join(contentDir, 'notes.json');
const projectsJsonPath = path.join(contentDir, 'projects.json');
const topicsJsonPath = path.join(contentDir, 'topics.json');
const categoriesYamlPath = path.join(contentSrcDir, 'categories.yaml');
const categoriesJsonPath = path.join(contentDir, 'categories.json');


// -------------------------- 辅助函数 --------------------------
/**
 * 确保目录存在（不存在则创建）
 * @param {string} filePath - 文件路径（用于推导目录）
 */
function ensureDirectoryExistence(filePath) {
  const dir = path.dirname(filePath);
  if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
}

/**
 * 安全处理数组（非数组则返回空数组）
 * @param {any} x - 输入值
 * @returns {any[]} 处理后的数组
 */
function safeArray(x) {
  return Array.isArray(x) ? x : [];
}

/**
 * 读取JSON数组文件（失败/不存在则返回默认值）
 * @param {string} path - 文件路径
 * @param {any[]} fallback - 默认值
 * @returns {any[]} 读取的数组
 */
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

/**
 * 读取YAML文件（失败/不存在则返回默认值）
 * @param {string} path - 文件路径
 * @param {object} fallback - 默认值
 * @returns {object} 读取的对象
 */
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

/**
 * 从文章URL提取子分类Key（用于匹配预定义分类）
 * 示例：/article/notes/Omics/genomics/bwa → 提取 "genomics"
 * @param {string} url - 文章URL
 * @returns {string} 子分类Key
 */
function getSubCategoryKeyFromUrl(url) {
  if (typeof url !== 'string') return '';
  const parts = url.replace(/^\/+/, '').split('/'); // 去除开头斜杠并拆分
  const articleIndex = parts[0] === 'article' ? 1 : 0; // 定位 "notes/projects/topics" 所在索引
  return parts[articleIndex + 2] || ''; // 子分类在 URL 中的位置（第4段）
}

/**
 * 生成文章的基础结构（标题+URL）
 * @param {object} article - 原始文章数据
 * @param {string} type - 类型（notes/projects/topics）
 * @returns {object} 格式化后的文章结构
 */
function formatArticle(article, type) {
  const relativePath = article?.relativePath || '';
  return {
    title: typeof article.title === 'string' ? article.title : relativePath.split('/').pop() || '未命名',
    articleUrl: `/article/${type}/${relativePath}`,
    wordCount: Number(article?.wordCount) || 0, // 保留字数用于分类统计
    date: article?.date || '', // 保留日期用于排序（可选）
    tags: safeArray(article?.tags).filter(t => typeof t === 'string') // 补充标签用于聚合
  };
}


// -------------------------- 数据归一化函数 --------------------------
/**
 * 归一化笔记的YAML配置（处理空值/类型错误）
 * @param {object} rawDef - 原始YAML配置
 * @returns {object} 归一化后的配置
 */
function normalizeNoteConfig(rawDef) {
  return {
    name: typeof rawDef?.name === 'string' ? rawDef.name : '',
    title: typeof rawDef?.title === 'string' ? rawDef.title : '',
    desc: typeof rawDef?.desc === 'string' ? rawDef.desc : '',
    date: typeof rawDef?.date === 'string' ? rawDef.date : '',
    categories: rawDef?.categories && typeof rawDef.categories === 'object' ? rawDef.categories : {}
  };
}

/**
 * 归一化项目/课题的YAML配置（处理空值/类型错误）
 * @param {object} rawDef - 原始YAML配置
 * @returns {object} 归一化后的配置
 */
function normalizeProjectTopicConfig(rawDef) {
  return {
    name: typeof rawDef?.name === 'string' ? rawDef.name : '',
    title: typeof rawDef?.title === 'string' ? rawDef.title : '',
    desc: typeof rawDef?.desc === 'string' ? rawDef.desc : '',
    date: typeof rawDef?.date === 'string' ? rawDef.date : '',
    tags: safeArray(rawDef?.tags).filter(t => typeof t === 'string'),
    github: typeof rawDef?.github === 'string' ? rawDef.github : '',
    doi: typeof rawDef?.doi === 'string' ? rawDef.doi : '',
    url: typeof rawDef?.url === 'string' ? rawDef.url : '',
    categories: rawDef?.categories && typeof rawDef.categories === 'object' ? rawDef.categories : {}
  };
}


// -------------------------- 核心分类构建函数 --------------------------
/**
 * 构建笔记的精细化分类（合并分类与文章）
 * @param {object[]} noteConfigs - 笔记的YAML配置列表
 * @param {object[]} noteArticles - 笔记的原始文章数据
 * @returns {object[]} 精细化分类后的笔记列表
 */
function buildDetailedNoteCategories(noteConfigs, noteArticles) {
  return noteConfigs.map(rawConfig => {
    const config = normalizeNoteConfig(rawConfig);
    if (!config.name) return null; // 过滤无名称的配置

    // 1. 筛选当前笔记分类下的所有文章（按 name 匹配）
    const currentArticles = noteArticles
      .filter(art => typeof art?.relativePath === 'string' && art.relativePath.split('/')[0] === config.name)
      .map(art => formatArticle(art, 'notes')); // 格式化文章结构

    // 2. 按预定义子分类拆分文章（生成精细化分类）
    const predefinedSubCats = Object.entries(config.categories); // [['genomics', '基因组学'], ...]
    const detailedSubCats = predefinedSubCats.map(([key, title]) => {
      // 筛选该子分类下的文章（通过URL匹配子分类Key）
      const catArticles = currentArticles.filter(art => getSubCategoryKeyFromUrl(art.articleUrl) === key);
      // 子分类级统计
      return {
        key, // 子分类标识（用于匹配）
        title, // 子分类显示名
        articles: catArticles.map(({ wordCount, date, ...rest }) => rest), // 移除统计字段，保留展示用数据
        stats: {
          postsCount: catArticles.length,
          totalWords: catArticles.reduce((sum, art) => sum + art.wordCount, 0),
          latestDate: catArticles.length 
            ? catArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date 
            : ''
        }
      };
    });

    // 3. 处理未被分类的文章（添加"未分类"子分类）
    const categorizedUrls = detailedSubCats.flatMap(cat => cat.articles.map(art => art.articleUrl));
    const uncategorizedArticles = currentArticles.filter(art => !categorizedUrls.includes(art.articleUrl));
    if (uncategorizedArticles.length > 0) {
      detailedSubCats.push({
        key: 'uncategorized',
        title: '未分类',
        articles: uncategorizedArticles.map(({ wordCount, date, ...rest }) => rest),
        stats: {
          postsCount: uncategorizedArticles.length,
          totalWords: uncategorizedArticles.reduce((sum, art) => sum + art.wordCount, 0),
          latestDate: uncategorizedArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date
        }
      });
    }

    // 4. 全局统计（当前笔记分类的总数据）
    const globalStats = {
      postsCount: currentArticles.length,
      totalWords: currentArticles.reduce((sum, art) => sum + art.wordCount, 0),
      latestDate: config.date || (currentArticles.length 
        ? currentArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date 
        : '')
    };

    // 5. 全局标签（所有文章标签去重）
    const globalTags = Array.from(
      new Set(currentArticles.flatMap(art => safeArray(art?.tags || []).filter(t => typeof t === 'string')))
    );

    // 6. 计算Root URL（按分类顺序取第一个有文章的子分类的第一篇文章）
    const rootUrl = detailedSubCats
      .find(cat => cat.stats.postsCount > 0)
      ?.articles[0]?.articleUrl || '';

    return {
      name: config.name,
      title: config.title || config.name,
      desc: config.desc,
      tags: globalTags,
      stats: globalStats,
      categories: detailedSubCats, // 核心：合并后的精细化分类
      root: rootUrl
    };
  }).filter(Boolean); // 过滤无效配置
}

/**
 * 构建项目/课题的精细化分类（合并分类与文章）
 * @param {object[]} ptConfigs - 项目/课题的YAML配置列表
 * @param {object[]} ptArticles - 项目/课题的原始文章数据
 * @param {string} type - 类型（project/topic）
 * @returns {object[]} 精细化分类后的项目/课题列表
 */
function buildDetailedProjectTopicCategories(ptConfigs, ptArticles, type) {
  return ptConfigs.map(rawConfig => {
    const config = normalizeProjectTopicConfig(rawConfig);
    if (!config.name) return null; // 过滤无名称的配置

    // 1. 筛选当前项目/课题下的所有文章（按名称归一化匹配）
    const nameKey = config.name.toLowerCase().replace(/\s+/g, ''); // 归一化名称（去空格转小写）
    const nameAlpha = config.name.toLowerCase().replace(/[^a-z0-9]+/g, ''); // 仅保留字母数字
    const currentArticles = ptArticles
      .filter(art => {
        if (typeof art?.relativePath !== 'string') return false;
        const artName = art.relativePath.split('/')[0].toLowerCase().replace(/\s+/g, '');
        const artAlpha = art.relativePath.split('/')[0].toLowerCase().replace(/[^a-z0-9]+/g, '');
        return [nameKey, nameAlpha].includes(artName) || [nameKey, nameAlpha].includes(artAlpha);
      })
      .map(art => formatArticle(art, `${type}s`)); // 格式化文章结构（URL用复数：projects/topics）

    // 2. 按预定义子分类拆分文章（生成精细化分类）
    const predefinedSubCats = Object.entries(config.categories);
    const detailedSubCats = predefinedSubCats.map(([key, title]) => {
      const catArticles = currentArticles.filter(art => getSubCategoryKeyFromUrl(art.articleUrl) === key);
      return {
        key,
        title,
        articles: catArticles.map(({ wordCount, date, ...rest }) => rest),
        stats: {
          postsCount: catArticles.length,
          totalWords: catArticles.reduce((sum, art) => sum + art.wordCount, 0),
          latestDate: catArticles.length 
            ? catArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date 
            : ''
        }
      };
    });

    // 3. 处理未被分类的文章
    const categorizedUrls = detailedSubCats.flatMap(cat => cat.articles.map(art => art.articleUrl));
    const uncategorizedArticles = currentArticles.filter(art => !categorizedUrls.includes(art.articleUrl));
    if (uncategorizedArticles.length > 0) {
      detailedSubCats.push({
        key: 'uncategorized',
        title: '未分类',
        articles: uncategorizedArticles.map(({ wordCount, date, ...rest }) => rest),
        stats: {
          postsCount: uncategorizedArticles.length,
          totalWords: uncategorizedArticles.reduce((sum, art) => sum + art.wordCount, 0),
          latestDate: uncategorizedArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date
        }
      });
    }

    // 4. 全局统计
    const globalStats = {
      postsCount: currentArticles.length,
      totalWords: currentArticles.reduce((sum, art) => sum + art.wordCount, 0),
      latestDate: config.date || (currentArticles.length 
        ? currentArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date 
        : '')
    };

    // 5. 计算Root URL（优先用配置的url，否则取第一个有文章的子分类）
    let rootUrl = '';
    if (config.url && config.url.startsWith('/article/')) {
      rootUrl = config.url;
    } else {
      rootUrl = detailedSubCats.find(cat => cat.stats.postsCount > 0)?.articles[0]?.articleUrl || '';
    }

    // 6. 附加项目/课题特有字段（GitHub/DOI）
    const extraFields = type === 'project' 
      ? { github: config.github } 
      : { doi: config.doi };

    return {
      name: config.name,
      title: config.title || config.name,
      desc: config.desc,
      tags: config.tags,
      stats: globalStats,
      categories: detailedSubCats, // 核心：合并后的精细化分类
      root: rootUrl,
      ...extraFields
    };
  }).filter(Boolean); // 过滤无效配置
}


// -------------------------- 主流程函数 --------------------------
/**
 * 生成精细化分类的JSON文件
 */
function generateCategoriesJson() {
  try {
    // 1. 读取原始数据
    const noteArticles = readJsonArray(notesJsonPath);
    const projectArticles = readJsonArray(projectsJsonPath);
    const topicArticles = readJsonArray(topicsJsonPath);
    const yamlConfig = readYaml(categoriesYamlPath);
    const { notes: noteConfigs, projects: projectConfigs, topics: topicConfigs } = {
      notes: safeArray(yamlConfig?.notes),
      projects: safeArray(yamlConfig?.projects),
      topics: safeArray(yamlConfig?.topics)
    };

    // 2. 构建各类型的精细化分类
    const detailedNotes = buildDetailedNoteCategories(noteConfigs, noteArticles);
    const detailedProjects = buildDetailedProjectTopicCategories(projectConfigs, projectArticles, 'project');
    const detailedTopics = buildDetailedProjectTopicCategories(topicConfigs, topicArticles, 'topic');

    // 3. 合并为最终结构（保持原有的"笔记/项目/课题"三级分类）
    const finalStructure = [
      {
        title: '笔记',
        items: detailedNotes.sort((a, b) => a.name.localeCompare(b.name, 'zh-Hans-CN')) // 按名称排序
      },
      {
        title: '项目',
        items: detailedProjects.sort((a, b) => a.name.localeCompare(b.name, 'zh-Hans-CN'))
      },
      {
        title: '课题',
        items: detailedTopics.sort((a, b) => a.name.localeCompare(b.name, 'zh-Hans-CN'))
      }
    ];

    // 4. 写入文件
    ensureDirectoryExistence(categoriesJsonPath);
    fs.writeFileSync(categoriesJsonPath, JSON.stringify(finalStructure, null, 2), 'utf-8');
    console.log(`Successfully generated: ${categoriesJsonPath}`);
  } catch (error) {
    console.error(`Failed to generate ${categoriesJsonPath}:`, error);
    process.exitCode = 1;
  }
}

/**
 * 主函数（脚本入口）
 */
function main() {
  console.log('Starting categories.json generation script...');
  generateCategoriesJson();
  console.log('categories.json generation complete.');
}

// 脚本直接运行时执行主函数
if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main();
}

// 导出函数供外部调用（可选）
export { generateCategoriesJson };