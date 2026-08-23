/**
 * 翻译器配置文件
 */

import path from 'path'

import { contentSrcDir, EN_SUFFIX } from '../../dataConfig.ts'

export const TRANSLATION_CONFIG = {
  SUPPORTED_EXTENSIONS: ['.md', '.yaml', '.yml'],

  /** 输出后缀（-en 是内容身份，与 dataConfig.EN_SUFFIX 单一来源） */
  OUTPUT_SUFFIX: EN_SUFFIX,

  EXCLUDE_PATTERNS: [`*${EN_SUFFIX}.*`, '*.bak', '*.tmp'],

  DEFAULT_OPTIONS: {
    force: false,
    skipExisting: true,
    recursive: true,
    filePatterns: ['*.md', '*.yaml', '*.yml'],
  },
}

/** 翻译目标：输入一律为 content-src 手写中文源，输出镜像写入 cache/en */
export const TRANSLATION_TARGETS = {
  CONTENT_SRC: contentSrcDir,

  CONFIG_FILES: {
    ABOUT: path.join(contentSrcDir, 'about.yaml'),
    CATEGORIES: path.join(contentSrcDir, 'categories.yaml'),
    RESOURCES: path.join(contentSrcDir, 'resources.yaml'),
  },
}
