/**
 * 脚本/工具的统一日志出口：所有终端输出走这里，保证前缀、级别与语言一致。
 * [INFO] 进度 / [OK] 成功结果 / [WARN] 可继续的告警（stderr）/
 * [ERROR] 失败（stderr）/ [FIX] 自动修复动作；== 标题 == 为阶段横幅。
 * 数据正文（列表/报表体）不带标签，仅状态与结论行带。
 */

type Stream = 'log' | 'warn' | 'error'
type Tag = 'INFO' | 'OK' | 'WARN' | 'ERROR' | 'FIX'

const STREAM_BY_TAG: Record<Tag, Stream> = {
  INFO: 'log',
  OK: 'log',
  WARN: 'warn',
  ERROR: 'error',
  FIX: 'log',
}

function emit(tag: Tag, message: string, detail?: unknown): void {
  const stream = STREAM_BY_TAG[tag]
  if (detail === undefined) console[stream](`[${tag}] ${message}`)
  else console[stream](`[${tag}] ${message}`, detail)
}

export function logInfo(message: string): void {
  emit('INFO', message)
}

export function logOk(message: string): void {
  emit('OK', message)
}

export function logWarn(message: string, detail?: unknown): void {
  emit('WARN', message, detail)
}

export function logError(message: string, detail?: unknown): void {
  emit('ERROR', message, detail)
}

export function logFix(message: string): void {
  emit('FIX', message)
}

export function logSection(title: string): void {
  console.log(`== ${title} ==`)
}
