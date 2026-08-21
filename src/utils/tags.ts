import type { Post } from '@/types'

export type TagLevel = 'lg' | 'md' | 'sm'

export interface TagCloudItem {
  name: string
  count: number
  level: TagLevel
}

/** 标签云：统计词频并按词频降序三等分档（lg/md/sm），输出保持 name 字母序 */
export function buildTagCloud(posts: Post[], levelCount = 3): TagCloudItem[] {
  const map = new Map<string, number>()
  posts.forEach((p) => p.tags.forEach((t) => map.set(t, (map.get(t) || 0) + 1)))
  const list = Array.from(map.entries())
    .map(([name, count]) => ({ name, count }))
    .sort((a, b) => a.name.localeCompare(b.name))
  const total = list.length
  const third = Math.ceil(total / levelCount)
  const rank = new Map(
    [...list].sort((a, b) => b.count - a.count).map((item, idx) => [item.name, idx]),
  )
  return list.map((item) => {
    const idx = rank.get(item.name) ?? 0
    return { ...item, level: idx < third ? 'lg' : idx < third * 2 ? 'md' : 'sm' }
  })
}
