/**
 * 分页窗口计算：返回含省略号占位的可见页码序列（1 与 total 恒为端点）。
 * total <= maxShow 时全量显示；窗口贴首/贴尾时自动收拢并补省略号。
 */
export function getVisiblePages(
  current: number,
  total: number,
  maxShow = 5,
): Array<number | '...'> {
  if (total <= maxShow) {
    return Array.from({ length: total }, (_, i) => i + 1)
  }

  const pages: number[] = [current]

  let i = 1
  while (pages.length < maxShow) {
    if (current - i >= 1 && !pages.includes(current - i)) {
      pages.push(current - i)
    }
    if (pages.length >= maxShow) break

    if (current + i <= total && !pages.includes(current + i)) {
      pages.push(current + i)
    }
    i++
  }

  if (!pages.includes(1)) {
    pages.pop()
    pages.push(1)
  }
  if (pages.length < maxShow && !pages.includes(total)) {
    pages.push(total)
  } else if (pages.length >= maxShow && !pages.includes(total)) {
    pages.pop()
    pages.push(total)
  }

  const sorted = pages.sort((a, b) => a - b)
  const result: Array<number | '...'> = []
  for (let idx = 0; idx < sorted.length; idx++) {
    if (idx > 0 && sorted[idx] - sorted[idx - 1] > 1) {
      result.push('...')
    }
    result.push(sorted[idx])
  }
  /** 首部被窗口挤出时用省略号占位（尾部恒有 total，无需补位） */
  if (result[0] !== 1) {
    result.unshift('...')
  }
  return result
}
