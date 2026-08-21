/**
 * 语言包结构基准：以 zh-CN 为基准，en-US 必须 `satisfies AppMessages`
 * （键缺失 / 类型不匹配在编译期报错，防止双语键漂移）。
 */
import zhCN from './locales/zh-CN'

export type AppMessages = typeof zhCN
