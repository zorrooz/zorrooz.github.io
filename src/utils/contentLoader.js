/**
 * 统一的内容加载工具函数
 * 根据当前语言动态加载对应的JSON或MD文件
 */

/**
 * 获取当前语言环境
 * @returns {string} 当前语言代码 ('zh-CN' 或 'en-US')
 */
export const getCurrentLocale = () => {
  return localStorage.getItem('locale') || 'zh-CN';
};

/**
 * 根据当前语言生成对应的文件名
 * @param {string} baseName - 基础文件名（不含后缀）
 * @param {string} extension - 文件扩展名（如 '.json', '.md'）
 * @returns {string} 对应的语言文件名
 */
export const getLocalizedFileName = (baseName, extension = '.json') => {
  const locale = getCurrentLocale();
  const isEnglish = locale === 'en-US';
  
  // 如果已经是英文文件，直接返回
  if (baseName.endsWith('-en')) {
    return `${baseName}${extension}`;
  }
  
  // 如果是中文文件，根据语言决定是否添加-en后缀
  return isEnglish ? `${baseName}-en${extension}` : `${baseName}${extension}`;
};

/**
 * 动态加载JSON内容文件
 * @param {string} fileName - 基础文件名（不含后缀）
 * @param {string} directory - 文件所在目录（相对于src/content）
 * @returns {Promise<Object>} 加载的JSON数据
 */
export const loadJsonContent = async (fileName, directory = '') => {
  const localizedFileName = getLocalizedFileName(fileName, '.json');
  
  try {
    // 使用Vite支持的import.meta.glob进行动态导入
    const modules = import.meta.glob('../content/**/*.json', { eager: true });
    
    // 构建完整的文件路径
    const fullPath = directory ? `../content/${directory}/${localizedFileName}` : `../content/${localizedFileName}`;
    
    // 查找匹配的模块
    const matchedKey = Object.keys(modules).find(key => key.includes(localizedFileName));
    
    if (matchedKey) {
      return modules[matchedKey].default || {};
    }
    
    throw new Error(`JSON file not found: ${localizedFileName}`);
  } catch (error) {
    console.error(`Failed to load JSON content: ${localizedFileName}`, error);
    
    // 回退到默认中文文件
    if (getCurrentLocale() === 'en-US') {
      try {
        const modules = import.meta.glob('../content/**/*.json', { eager: true });
        const fallbackKey = Object.keys(modules).find(key => key.includes(`${fileName}.json`));
        
        if (fallbackKey) {
          return modules[fallbackKey].default || {};
        }
      } catch (fallbackError) {
        console.error(`Failed to load fallback JSON content: ${fileName}.json`, fallbackError);
      }
    }
    
    return {};
  }
};

/**
 * 动态加载MD文件内容
 * @param {string} filePath - MD文件路径（相对于src/content-src）
 * @returns {Promise<string>} MD文件内容
 */
export const loadMarkdownContent = async (filePath) => {
  const locale = getCurrentLocale();
  const isEnglish = locale === 'en-US';
  
  // 构建所有可能的文件路径组合（按优先级排序）
  const possiblePaths = [];
  
  // 1. 当前语言对应的文件（最高优先级）
  const localizedPath = isEnglish ? filePath.replace('.md', '-en.md') : filePath;
  possiblePaths.push(localizedPath);
  
  // 2. 原始文件（中等优先级）
  possiblePaths.push(filePath);
  
  // 3. 跨语言回退（最低优先级）
  if (isEnglish) {
    // 英文环境回退到中文
    const chinesePath = filePath.replace('-en.md', '.md');
    if (chinesePath !== filePath && chinesePath !== localizedPath) {
      possiblePaths.push(chinesePath);
    }
  } else {
    // 中文环境回退到英文
    const englishPath = filePath.replace('.md', '-en.md');
    if (englishPath !== filePath && englishPath !== localizedPath) {
      possiblePaths.push(englishPath);
    }
  }
  
  try {
    // 使用Vite的glob导入模式
    const markdownModules = import.meta.glob('../content-src/**/*.md', { 
      query: '?raw', 
      import: 'default', 
      eager: false 
    });
    
    // 按优先级顺序尝试所有可能的路径
    for (const path of possiblePaths) {
      const searchPaths = [
        `../content-src/${path}`,
        `../content-src/${path}?raw`
      ];
      
      const matchedKey = Object.keys(markdownModules).find(key => 
        searchPaths.some(searchPath => key.endsWith(searchPath))
      );
      
      if (matchedKey) {
        const content = await markdownModules[matchedKey]();
        return content;
      }
    }
    
    throw new Error(`Markdown file not found for any of: ${possiblePaths.join(', ')}`);
  } catch (error) {
    console.error(`Failed to load markdown content: ${localizedPath}`, error);
    
    // 回退到原始文件
    if (isEnglish) {
      try {
        const markdownModules = import.meta.glob('../content-src/**/*.md', { 
          query: '?raw', 
          import: 'default', 
          eager: false 
        });
        
        const fallbackPaths = [
          `../content-src/${filePath}`,
          `../content-src/${filePath}?raw`
        ];
        
        const fallbackKey = Object.keys(markdownModules).find(key => 
          fallbackPaths.some(path => key.endsWith(path))
        );
        
        if (fallbackKey) {
          const content = await markdownModules[fallbackKey]();
          return content;
        }
      } catch (fallbackError) {
        console.error(`Failed to load fallback markdown content: ${filePath}`, fallbackError);
      }
    }
    
    return '# Content Not Available\n\nThe requested content could not be loaded.';
  }
};

/**
 * 批量加载多个JSON文件
 * @param {Array<string>} fileNames - 文件名数组
 * @param {string} directory - 文件所在目录
 * @returns {Promise<Object>} 合并的数据对象
 */
export const loadMultipleJsonContents = async (fileNames, directory = '') => {
  const results = {};
  
  for (const fileName of fileNames) {
    try {
      const data = await loadJsonContent(fileName, directory);
      results[fileName] = data;
    } catch (error) {
      console.error(`Failed to load ${fileName}`, error);
      results[fileName] = {};
    }
  }
  
  return results;
};

// 常用内容文件的快捷加载方法
export const loadCategories = () => loadJsonContent('categories');
export const loadPosts = () => loadJsonContent('posts');
export const loadNotes = () => loadJsonContent('notes');
export const loadProjects = () => loadJsonContent('projects');
export const loadTopics = () => loadJsonContent('topics');
export const loadTags = () => loadJsonContent('tags');
export const loadAbout = () => loadJsonContent('about');
export const loadResources = () => loadJsonContent('resources');

export default {
  getCurrentLocale,
  getLocalizedFileName,
  loadJsonContent,
  loadMarkdownContent,
  loadMultipleJsonContents,
  loadCategories,
  loadPosts,
  loadNotes,
  loadProjects,
  loadTopics,
  loadTags,
  loadAbout,
  loadResources
};