/**
 * 翻译器配置文件
 */

// 支持的翻译配置
export const TRANSLATION_CONFIG = {
  // 支持的文件类型
  SUPPORTED_EXTENSIONS: ['.md', '.yaml', '.yml'],
  
  // 输出文件后缀
  OUTPUT_SUFFIX: '-en',
  
  // 排除的文件模式
  EXCLUDE_PATTERNS: ['*-en.*', '*.bak', '*.tmp'],
  
  // 默认翻译选项
  DEFAULT_OPTIONS: {
    force: false,
    skipExisting: true,
    recursive: true,
    filePatterns: ['*.md', '*.yaml', '*.yml']
  },
  
  // 翻译策略
  STRATEGIES: {
    INCREMENTAL: 'incremental', // 增量翻译（默认）
    FORCE: 'force',            // 强制重新翻译
    NEW_ONLY: 'new_only'       // 仅翻译新文件
  }
};

// 翻译目标配置
export const TRANSLATION_TARGETS = {
  // 主要翻译目录
  CONTENT_SRC: 'src/content-src',
  
  // 子目录配置
  SUBDIRECTORIES: {
    NOTES: 'src/content-src/notes',
    PROJECTS: 'src/content-src/projects', 
    TOPICS: 'src/content-src/topics',
    TEST: 'src/content-src/test'
  },
  
  // 配置文件
  CONFIG_FILES: {
    ABOUT: 'src/content-src/about.yaml',
    CATEGORIES: 'src/content-src/categories.yaml',
    RESOURCES: 'src/content-src/resources.yaml'
  }
};

// 日志配置
export const LOG_CONFIG = {
  LEVELS: {
    ERROR: 'error',
    WARN: 'warn', 
    INFO: 'info',
    DEBUG: 'debug'
  },
  
  // 日志前缀
  PREFIXES: {
    SUCCESS: '[INFO]',
    WARNING: '[WARN]',
    ERROR: '[ERROR]',
    PROCESSING: '[INFO]',
    SKIP: '[INFO]'
  }
};