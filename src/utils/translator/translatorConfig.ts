/**
 * 翻译器配置文件
 */

import path from 'path'

import { contentSrcDir } from '../dataConfig.ts'

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

// 翻译目标配置（数据目录由 dataConfig 统一提供）。
// 输入一律是 content-src 的手写中文源，输出由 translator 镜像写入 cache/en（第二层）。
export const TRANSLATION_TARGETS = {
  // 主要翻译目录（notes 文章 + yaml 配置文件）
  CONTENT_SRC: contentSrcDir,

  // 配置文件（手写中文源）
  CONFIG_FILES: {
    ABOUT: path.join(contentSrcDir, 'about.yaml'),
    CATEGORIES: path.join(contentSrcDir, 'categories.yaml'),
    RESOURCES: path.join(contentSrcDir, 'resources.yaml'),
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
