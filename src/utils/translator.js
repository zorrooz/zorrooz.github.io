import OpenAI from "openai";
import fs from "fs/promises";
import path from "path";
import apiKeys from "../config/apiKeys.js";

const openai = new OpenAI({
  baseURL: 'https://api.deepseek.com',
  apiKey: apiKeys.deepseek
});

/**
 * 翻译文本内容
 * @param {string} text - 要翻译的文本
 * @param {string} fileType - 文件类型 ('md' 或 'yaml')
 * @returns {Promise<string>} 翻译后的文本
 */
async function translateText(text, fileType = 'md') {
  try {
    const systemPrompt = fileType === 'yaml' 
      ? "你是一个专业翻译器。请将以下YAML文件内容翻译为英文，要求：\n1. 严格保持YAML格式结构（键名不翻译，只翻译值）\n2. 保持缩进、标点等格式不变\n3. 不添加任何解释、注释或额外内容\n4. 只输出翻译结果"
      : "你是一个专业翻译器。请将以下Markdown内容翻译为英文，要求：\n1. 严格保持原始格式（包括Markdown语法、代码块、换行、缩进等）\n2. 不添加任何解释、注释或额外内容\n3. 只输出翻译结果";

    const completion = await openai.chat.completions.create({
      messages: [
        { role: "system", content: systemPrompt },
        { role: "user", content: text }
      ],
      model: "deepseek-chat",
      temperature: 0.3
    });

    return completion.choices[0].message.content.trim();
  } catch (error) {
    console.error("翻译时发生错误:", error.message);
    throw error;
  }
}

/**
 * 检查文件是否需要翻译（增量翻译）
 * @param {string} sourcePath - 源文件路径
 * @param {string} targetPath - 目标文件路径
 * @returns {Promise<boolean>} 是否需要翻译
 */
async function needsTranslation(sourcePath, targetPath) {
  try {
    // 检查目标文件是否存在
    await fs.access(targetPath);
    
    // 获取两个文件的修改时间
    const sourceStats = await fs.stat(sourcePath);
    const targetStats = await fs.stat(targetPath);
    
    // 如果源文件比目标文件新，则需要重新翻译
    return sourceStats.mtime > targetStats.mtime;
  } catch (error) {
    // 目标文件不存在，需要翻译
    return true;
  }
}

/**
 * 翻译单个文件
 * @param {string} inputFilePath - 输入文件路径
 * @param {Object} options - 选项
 * @returns {Promise<void>}
 */
async function translateFile(inputFilePath, options = {}) {
  const {
    force = false, // 强制重新翻译
    skipExisting = true, // 跳过已存在的翻译文件
    outputSuffix = '-en' // 输出文件后缀
  } = options;

  try {
    // 检查文件扩展名
    const ext = path.extname(inputFilePath).toLowerCase();
    if (!['.md', '.yaml', '.yml'].includes(ext)) {
      console.log(`[WARN] 跳过不支持的文件类型: ${inputFilePath}`);
      return;
    }

    // 生成输出路径
    const basename = path.basename(inputFilePath, ext);
    const outputPath = path.join(
      path.dirname(inputFilePath),
      `${basename}${outputSuffix}${ext}`
    );

    // 检查是否需要翻译
    if (skipExisting && !force) {
      const shouldTranslate = await needsTranslation(inputFilePath, outputPath);
      if (!shouldTranslate) {
        console.log(`[INFO] 跳过已翻译文件: ${inputFilePath}`);
        return null; // 返回null而不是undefined
      }
    }

    // 读取文件内容
    const content = await fs.readFile(inputFilePath, "utf-8");
    
    // 翻译
    console.log(`[INFO] 正在翻译文件: ${inputFilePath}`);
    const fileType = ext === '.yaml' || ext === '.yml' ? 'yaml' : 'md';
    const translated = await translateText(content, fileType);

    // 写入翻译后文件
    await fs.writeFile(outputPath, translated, "utf-8");
    console.log(`[INFO] 翻译完成: ${outputPath}`);
    
    return outputPath;
  } catch (error) {
    console.error(`[ERROR] 处理文件时出错 (${inputFilePath}):`, error.message);
    throw error;
  }
}

/**
 * 批量翻译目录中的文件
 * @param {string} directoryPath - 目录路径
 * @param {Object} options - 选项
 * @returns {Promise<Array>} 翻译的文件列表
 */
async function translateDirectory(directoryPath, options = {}) {
  const {
    recursive = true, // 是否递归搜索子目录
    filePatterns = ['*.md', '*.yaml', '*.yml'], // 文件模式
    excludePatterns = ['*-en.*'], // 排除模式
    ...translateOptions
  } = options;

  try {
    const files = [];
    
    // 递归搜索文件
    async function searchFiles(dir) {
      const entries = await fs.readdir(dir, { withFileTypes: true });
      
      for (const entry of entries) {
        const fullPath = path.join(dir, entry.name);
        
        if (entry.isDirectory() && recursive) {
          await searchFiles(fullPath);
        } else if (entry.isFile()) {
          // 检查文件是否匹配模式且不被排除
          const matchesPattern = filePatterns.some(pattern => {
            const regex = new RegExp(pattern.replace('*', '.*').replace('.', '\\.'));
            return regex.test(entry.name);
          });
          
          const isExcluded = excludePatterns.some(pattern => {
            const regex = new RegExp(pattern.replace('*', '.*').replace('.', '\\.'));
            return regex.test(entry.name);
          });
          
          if (matchesPattern && !isExcluded) {
            files.push(fullPath);
          }
        }
      }
    }
    
    await searchFiles(directoryPath);
    
    console.log(`[INFO] 找到 ${files.length} 个需要翻译的文件`);
    
    // 批量翻译文件
    const results = [];
    for (const file of files) {
      try {
        const result = await translateFile(file, translateOptions);
        if (result) results.push(result);
      } catch (error) {
        console.error(`[ERROR] 翻译失败: ${file}`, error.message);
      }
    }
    
    console.log(`[INFO] 批量翻译完成，成功翻译 ${results.length} 个文件`);
    return results;
  } catch (error) {
    console.error('[ERROR] 批量翻译时出错:', error.message);
    throw error;
  }
}

/**
 * 翻译管理器 - 支持多种翻译策略
 */
class TranslationManager {
  constructor(options = {}) {
    this.options = options;
  }
  
  /**
   * 翻译单个文件或目录
   * @param {string} targetPath - 目标路径
   * @param {Object} options - 选项
   */
  async translate(targetPath, options = {}) {
    const mergedOptions = { ...this.options, ...options };
    
    try {
      const stats = await fs.stat(targetPath);
      
      if (stats.isDirectory()) {
        return await translateDirectory(targetPath, mergedOptions);
      } else {
        return await translateFile(targetPath, mergedOptions);
      }
    } catch (error) {
      console.error("[ERROR] 路径检查失败:", error.message);
      throw error;
    }
  }
  
  /**
   * 强制重新翻译所有文件
   * @param {string} targetPath - 目标路径
   */
  async forceTranslate(targetPath) {
    return this.translate(targetPath, { force: true, skipExisting: false });
  }
  
  /**
   * 仅翻译新文件（跳过所有已存在的翻译）
   * @param {string} targetPath - 目标路径
   */
  async translateNewOnly(targetPath) {
    return this.translate(targetPath, { skipExisting: true });
  }
}

// 导出函数和类
export { 
  translateText, 
  translateFile, 
  translateDirectory, 
  TranslationManager,
  needsTranslation 
};

// 默认导出管理器实例
export default new TranslationManager();

// CLI 使用示例
async function main() {
  const manager = new TranslationManager();
  
  // 示例用法
  try {
    // 翻译整个 content-src 目录（增量模式）
    console.log('[INFO] 开始增量翻译 content-src 目录...');
    const results = await manager.translate('../content-src');
    console.log(`[INFO] 翻译完成，处理了 ${results.length} 个文件`);
    
    // 强制重新翻译所有文件
    // await manager.forceTranslate('../content-src');
    
    // 仅翻译新文件
    // await manager.translateNewOnly('../content-src');
    
  } catch (error) {
    console.error('[ERROR] 翻译过程出错:', error.message);
  }
}

// 如果直接运行此文件，执行main函数
if (import.meta.url === `file://${process.argv[1]}`) {
  main();
}