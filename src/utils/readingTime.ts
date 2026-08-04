/**
 * 阅读时间估算（全站唯一实现，首页列表与文章页共用）。
 * 统一按 300 词/分钟估算（中文约 300 字/分钟）。
 */
export function readingTimeMinutes(wordCount: number): number {
  return Math.max(1, Math.round(Number(wordCount) / 300))
}
