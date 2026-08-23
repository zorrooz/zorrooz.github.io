import type { Post } from '@/types'

export interface TagCloudItem {
  name: string
  count: number
}

/** 标签云：统计词频，输出保持 name 字母序 */
export function buildTagCloud(posts: Post[]): TagCloudItem[] {
  const map = new Map<string, number>()
  posts.forEach((p) => p.tags.forEach((t) => map.set(t, (map.get(t) || 0) + 1)))
  return Array.from(map.entries())
    .map(([name, count]) => ({ name, count }))
    .sort((a, b) => a.name.localeCompare(b.name))
}
