/**
 * 翻译器配置文件
 */

import path from 'path'

import { contentSrcDir, EN_SUFFIX } from '../dataConfig.ts'

// 支持的翻译配置
export const TRANSLATION_CONFIG = {
  // 支持的文件类型
  SUPPORTED_EXTENSIONS: ['.md', '.yaml', '.yml'],

  // 输出文件后缀（-en 是内容身份，与 dataConfig.EN_SUFFIX 单一来源）
  OUTPUT_SUFFIX: EN_SUFFIX,

  // 排除的文件模式
  EXCLUDE_PATTERNS: [`*${EN_SUFFIX}.*`, '*.bak', '*.tmp'],

  // 默认翻译选项
  DEFAULT_OPTIONS: {
    force: false,
    skipExisting: true,
    recursive: true,
    filePatterns: ['*.md', '*.yaml', '*.yml'],
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
