// scripts/lib/tagMapping.ts — tag-mapping.json 读取共享实现（translator 与 tagMerger 共用）
import fs from 'fs'

export interface TagMapping {
  locale: string
  generatedAt: string
  mapping: Record<string, string>
  translation?: Record<string, string>
}

/** 读取 tag-mapping.json（完整结构）；缺失或损坏返回 null */
export function loadTagMapping(mappingPath: string): TagMapping | null {
  try {
    if (!fs.existsSync(mappingPath)) return null
    const parsed = JSON.parse(fs.readFileSync(mappingPath, 'utf-8')) as TagMapping
    if (!parsed?.mapping || typeof parsed.mapping !== 'object') return null
    return parsed
  } catch {
    return null
  }
}
