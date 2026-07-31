/**
 * 翻译器配置文件
 */

import path from 'path'

// 项目根目录：translatorConfig.ts 位于 src/utils/translator/
const projectRoot = path.resolve(import.meta.dirname, '../../..')

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
    filePatterns: ['*.md', '*.yaml', '*.yml'],
  },

  // 翻译策略
  STRATEGIES: {
    INCREMENTAL: 'incremental', // 增量翻译（默认）
    FORCE: 'force', // 强制重新翻译
    NEW_ONLY: 'new_only', // 仅翻译新文件
  },
}

// 翻译目标配置（基于项目根的绝对路径，不依赖 cwd）
export const TRANSLATION_TARGETS = {
  // 主要翻译目录
  CONTENT_SRC: path.join(projectRoot, 'src/content-src'),

  // 子目录配置
  SUBDIRECTORIES: {
    NOTES: path.join(projectRoot, 'src/content-src/notes'),
    PROJECTS: path.join(projectRoot, 'src/content-src/projects'),
    TOPICS: path.join(projectRoot, 'src/content-src/topics'),
    TEST: path.join(projectRoot, 'src/content-src/test'),
  },

  // 配置文件
  CONFIG_FILES: {
    ABOUT: path.join(projectRoot, 'src/content-src/about.yaml'),
    CATEGORIES: path.join(projectRoot, 'src/content-src/categories.yaml'),
    RESOURCES: path.join(projectRoot, 'src/content-src/resources.yaml'),
  },
}

// 日志配置
export const LOG_CONFIG = {
  LEVELS: {
    ERROR: 'error',
    WARN: 'warn',
    INFO: 'info',
    DEBUG: 'debug',
  },

  // 日志前缀
  PREFIXES: {
    SUCCESS: '[INFO]',
    WARNING: '[WARN]',
    ERROR: '[ERROR]',
    PROCESSING: '[INFO]',
    SKIP: '[INFO]',
  },
}
