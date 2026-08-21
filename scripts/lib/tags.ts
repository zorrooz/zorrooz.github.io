// scripts/lib/tags.ts — 标签计数 / 排序共享实现（generateTags 与 tagMerger 共用）
/** 遍历元素：字符串 trim 后计数，空串与非字符串跳过 */
export function countTags(iterable: Iterable<unknown>): Map<string, number> {
  const map = new Map<string, number>()
  for (const raw of iterable) {
    const name = typeof raw === 'string' ? raw.trim() : ''
    if (!name) continue
    map.set(name, (map.get(name) || 0) + 1)
  }
  return map
}

/** 按标签名 zh-Hans-CN localeCompare 排序 */
export function sortTagsByName(
  map: Map<string, number>,
): Array<{ name: string; count: number }> {
  return Array.from(map.entries())
    .map(([name, count]) => ({ name, count }))
    .sort((a, b) => a.name.localeCompare(b.name, 'zh-Hans-CN'))
}
