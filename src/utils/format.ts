/** 数字紧凑格式化：>=100 万显示 M、>=1000 显示 K，非整时保留 1 位小数 */
export function formatCompactNumber(n: number): string {
  if (n >= 1_000_000) return (n / 1_000_000).toFixed(n % 1_000_000 ? 1 : 0) + 'M'
  if (n >= 1_000) return (n / 1_000).toFixed(n % 1_000 ? 1 : 0) + 'K'
  return String(n)
}
